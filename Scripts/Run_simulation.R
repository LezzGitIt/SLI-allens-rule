## Run the SLI methods-paper simulation and export all results needed by
## Allens_methods_sim.qmd and supplementary_info.qmd to Derived/Rds/simulation_results.rds.
##
## Run this script (source() or Rscript) any time the parameter grid, sample
## size, or evaluation logic changes. The downstream .qmd files read the .rds
## for every simulation-derived statistic they quote; they do still call
## gen_ex_data() to redraw the small illustrative example figures.

suppressPackageStartupMessages({
  library(tidyverse)
  library(smatr)
  library(MASS)
  library(janitor)
})

source("Scripts/Key_allometry_fns.R")

# Global settings -----------------------------------------------------------
sma_or_ma <- "SMA"
log_ratio <- TRUE   # TRUE = log(A/S) and log(A²/S); FALSE = log(A)/log(S) and log(A)²/log(S)
N_ind     <- "3,000"

r_12     <- c(.3, .45, .6)
r_13     <- c(0, -.15, -.3, -.5, -.7)
r_23     <- r_13
b_avg_12 <- c(0.22, 0.33, 0.44)

# Parameter matrix ------------------------------------------------------------
# Main grid: Longer / Fatter / Proportionally smaller -- three realizations of Bergmann's rule (temperature never increases mass or appendage length).
Shape <- c("Longer", "Proportionally smaller", "Fatter")
Shape <- setNames(Shape, Shape)

Parms_mat <- expand_grid(
  b_avg_12 = b_avg_12,
  r_12 = r_12,
  r_13 = r_13,
  r_23 = r_23
) %>%
  mutate(
    Temp_eff = case_when(
      r_13 > r_23 ~ Shape[1],
      r_13 == r_23 ~ Shape[2],
      r_13 < r_23 ~ Shape[3]
    ),
    Strength = abs(r_13 - r_23)
  ) %>%
  mutate(Scaling = case_when(
    b_avg_12 < 0.33 ~ "Hypoallometry",
    b_avg_12 > 0.33 ~ "Hyperallometry",
    near(b_avg_12, 0.33) ~ "Isometry"
  ))

# Proportionally larger: mirror of the negative "Proportionally smaller" diagonal -- temperature increases wing and mass equally (rather than decreasing both), i.e. species get uniformly bigger with warming with no true shape change. Simulated analogue of empirical "Inverse Bergmann's" species.
r_temp_pos <- c(.15, .3, .5, .7)

Parms_big <- expand_grid(
  b_avg_12 = b_avg_12,
  r_12     = r_12,
  r_temp   = r_temp_pos
) %>%
  mutate(r_13 = r_temp, r_23 = r_temp, r_temp = NULL,
         Temp_eff = "Proportionally larger",
         Strength = abs(r_13 - r_23),
         Scaling  = case_when(
           b_avg_12 < 0.33 ~ "Hypoallometry",
           b_avg_12 > 0.33 ~ "Hyperallometry",
           near(b_avg_12, 0.33) ~ "Isometry"
         ))

Parms_mat2 <- bind_rows(Parms_mat, Parms_big) %>%
  mutate(
    Temp_eff = factor(Temp_eff,
                      levels = c("Fatter", "Proportionally smaller", "Proportionally larger", "Longer")),
    Scaling = factor(Scaling,
                     levels = c("Hypoallometry", "Isometry", "Hyperallometry"))
  ) %>%
  relocate(Scaling, .before = b_avg_12)

Parms_mat3 <- Parms_mat2 %>% filter(!(r_13 == 0 & r_23 == 0))
N_spp <- nrow(Parms_mat3)

# Generate data ---------------------------------------------------------------
Cols <- Parms_mat3 %>% dplyr::select(starts_with(c("b_", "r_")))
N <- as.numeric(str_remove(N_ind, ","))

### gen_data() draws fresh individuals per species, and gen_cov_mat(vary_sd = TRUE) redraws each species' sd_log_morph from Uniform(0.05, 0.09), so this script produced different results on every run. Seed it so simulation_results.rds is reproducible from source.
# Regenerating the .rds will change every simulation-derived number the manuscript quotes, once. After that they are fixed.
set.seed(20260712)

df_morph_l <- pmap(Cols, \(...) gen_data(..., n = N, meas_error = 0,
                                         transient_error_mass = 0,
                                         transient_error_append = 0))

# Validate simulated direction via SMA intercepts ------------------------------
extract_sma_intercepts <- function(df) {
  mod_temp_bin     <- run_sma_mod(df, interaction = FALSE)
  mod_temp_bin_int <- run_sma_mod(df, interaction = TRUE)
  mod_parms     <- format_sma_parms(mod_temp_bin)
  mod_parms_int <- format_sma_parms(mod_temp_bin_int)
  tibble(
    cor_allometry     = cor(mod_parms$Temp_inc, mod_parms$elevation),
    cor_allometry_int = cor(mod_parms_int$Temp_inc, mod_parms_int$elevation)
  )
}

sma_intercepts  <- map(df_morph_l, extract_sma_intercepts) %>% list_rbind()
sma_intercepts2 <- bind_cols(Parms_mat3, sma_intercepts)

# Correlation between binned temperature and the SMA intercept, under both a
# fixed-slope and slope-varies-with-temperature model (@fig-cor-allometry-values
# in supplementary_info.qmd).
Parms_temp_bs <- sma_intercepts2 %>%
  filter(Scaling != "Inverse") %>%
  pivot_longer(cols = c(cor_allometry, cor_allometry_int),
               names_to = "SMA_mod", values_to = "Correlation") %>%
  mutate(SMA_mod = if_else(SMA_mod == "cor_allometry", "No interaction", "Interaction"))

# Species with the wrong-signed correlation are excluded from further analysis.
Sim_fail <- sma_intercepts2 %>%
  left_join(sma_intercepts2) %>%
  filter(
    Temp_eff == "Fatter"  & cor_allometry_int > -.2 |
    Temp_eff == "Fatter"  & cor_allometry     > -.2 |
    Temp_eff == "Longer"  & cor_allometry_int <  .2 |
    Temp_eff == "Longer"  & cor_allometry     <  .2
  )

# Fit all six methods to each simulated species --------------------------------
generate_metrics <- function(Sim_df) {
  sma_mod  <- sma(Append_log ~ Mass_log, data = Sim_df, method = sma_or_ma)
  est_b_sma <- coef(sma_mod)["slope"]
  Ols_mod  <- lm(Append_log ~ Mass_log, data = Sim_df, na.action = na.exclude)
  est_b_ols <- coef(Ols_mod)["Mass_log"]
  Sim_df <- Sim_df %>%
    mutate(resid_ols = residuals(Ols_mod),
           Append_mass  = if (log_ratio) Append_log - Mass_log    else Append_log / Mass_log,
           Append2_mass = if (log_ratio) 2*Append_log - Mass_log  else Append_log^2 / Mass_log) %>%
    calc_sli(b_sli = 0.33,    rename_col = "sli_isometry") %>%
    calc_sli(b_sli = est_b_sma, rename_col = "sli_estimated")
  Sim_df_s <- Sim_df %>% mutate(across(where(is.numeric), scale))
  list(Sim_df_s = Sim_df_s, coefs = tibble(est_b_sma, est_b_ols))
}

extract_coefs <- function(Sim_df_s, coefs) {
  ols_resid_app <- lm(resid_ols    ~ Temp_inc, data = Sim_df_s)
  ryding_app    <- lm(Append_log   ~ Mass_log + Temp_inc, data = Sim_df_s,
                      na.action = na.exclude)
  ratio_app     <- lm(Append_mass  ~ Temp_inc, data = Sim_df_s, na.action = na.exclude)
  ratio2_app    <- lm(Append2_mass ~ Temp_inc, data = Sim_df_s, na.action = na.exclude)
  sli_iso_app   <- lm(sli_isometry  ~ Temp_inc, data = Sim_df_s, na.action = na.exclude)
  sli_est_app   <- lm(sli_estimated ~ Temp_inc, data = Sim_df_s, na.action = na.exclude)
  tibble(
    coef_sli_iso_app   = coef(sli_iso_app)["Temp_inc"],
    coef_sli_est_app   = coef(sli_est_app)["Temp_inc"],
    coef_ols_resid_app = coef(ols_resid_app)["Temp_inc"],
    coef_ryding_app    = coef(ryding_app)["Temp_inc"],
    coef_ratio         = coef(ratio_app)["Temp_inc"],
    coef_ratio2        = coef(ratio2_app)["Temp_inc"],
    est_b_sma          = coefs$est_b_sma,
    est_b_ols          = coefs$est_b_ols,
    # Per-species Pearson correlation between each method's individual-level
    # metric and body mass (@fig-mass-cor-sim in supplementary_info.qmd),
    # confirming in simulated data the mechanical ratio/mass dependence shown
    # empirically in the main text (@tbl-ratio-mass-summary). Mass-as-covariate
    # is excluded: it has no individual-level shape metric to correlate, only
    # a model coefficient.
    cor_ratio          = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$Append_mass)),
    cor_ratio2         = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$Append2_mass)),
    cor_ols_resid      = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$resid_ols)),
    cor_sli_iso        = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$sli_isometry)),
    cor_sli_est        = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$sli_estimated))
  )
}

Parms_tbl <- map(df_morph_l, \(df) {
  m <- generate_metrics(df)
  extract_coefs(Sim_df_s = m$Sim_df_s, coefs = m$coefs)
}) %>% list_rbind()

Parms_tbl2 <- bind_cols(Parms_mat3, Parms_tbl)
Parms_tbl3 <- Parms_tbl2 %>% anti_join(Sim_fail)

x_labs <- c(
  "Sli_est"   = "SLI estimated",
  "Sli_iso"   = "SLI isometry",
  "Ols_resid" = "OLS residuals",
  "Ryding"    = "Mass as covariate",
  "Ratio2"    = "Appendage² / mass",
  "Ratio"     = "Appendage / mass"
)

Parms_tbl4 <- Parms_tbl3 %>%
  pivot_longer(
    cols = c(coef_sli_iso_app, coef_sli_est_app, coef_ols_resid_app,
             coef_ryding_app, coef_ratio, coef_ratio2),
    names_to  = "Model",
    values_to = "b_temp_inc"
  ) %>%
  mutate(
    Model = str_remove_all(Model, "coef_|_app"),
    Model = str_to_sentence(Model)
  )

# Per-species correlation between each method's metric and body mass ---------
Mass_cor_tbl <- Parms_tbl3 %>%
  pivot_longer(
    cols = c(cor_ratio, cor_ratio2, cor_ols_resid, cor_sli_iso, cor_sli_est),
    names_to  = "Model",
    values_to = "r_mass"
  ) %>%
  mutate(
    Model = str_remove_all(Model, "cor_"),
    Model = str_to_sentence(Model)
  )

# Evaluation --------------------------------------------------------------
Parms_tbl5 <- Parms_tbl4 %>%
  mutate(b_dir = if_else(b_temp_inc < 0, "Neg", "Pos")) %>%
  mutate(Est_correct = case_when(
    Temp_eff == "Longer"      & b_dir == "Pos" ~ TRUE,
    Temp_eff == "Fatter"      & b_dir == "Neg" ~ TRUE,
    Temp_eff == "Proportionally smaller" ~ NA,
    Temp_eff == "Proportionally larger"  ~ NA,
    .default = FALSE
  ), .by = Model)

Eval_tbl <- Parms_tbl5 %>%
  summarize(Prop_correct = sum(Est_correct, na.rm = TRUE) / n(),
            .by = c(Model, Temp_eff)) %>%
  mutate(Prop_correct = round(Prop_correct, 2)) %>%
  filter(!Temp_eff %in% c("Proportionally smaller", "Proportionally larger")) %>%
  arrange(Model)

Proportional_methods_eval <- Parms_tbl5 %>%
  filter(Temp_eff == "Proportionally smaller") %>%
  summarize(Prop_pos = sum(b_dir == "Pos") / n(), .by = c(Model, Temp_eff)) %>%
  mutate(Prop_neg = 1 - Prop_pos) %>%
  pivot_longer(cols = c(Prop_pos, Prop_neg),
               names_to = "Direction", values_to = "Proportion") %>%
  mutate(Direction = str_remove(Direction, "Prop_"),
         Model = str_replace_all(Model, x_labs))

Proportional_true_eff_tbl <- Parms_temp_bs %>%
  filter(Temp_eff == "Proportionally smaller") %>%
  mutate(True_temp_eff = if_else(Correlation > 0, "Longer", "Fatter")) %>%
  janitor::tabyl(True_temp_eff)

Prop_fatter <- Proportional_true_eff_tbl %>%
  filter(True_temp_eff == "Fatter") %>%
  mutate(percent = round(percent, 2)) %>%
  pull(percent)
Prop_longer <- 1 - Prop_fatter

# "Proportionally larger" scenario: same "no true shape signal" evaluation as Proportionally smaller, just mirrored to the positive-r_13/r_23 diagonal (species getting bigger, not smaller). Kept as a parallel/duplicate computation rather than generalizing Proportional_methods_eval, so those objects remain unaffected.
Bigger_methods_eval <- Parms_tbl5 %>%
  filter(Temp_eff == "Proportionally larger") %>%
  summarize(Prop_pos = sum(b_dir == "Pos") / n(), .by = c(Model, Temp_eff)) %>%
  mutate(Prop_neg = 1 - Prop_pos) %>%
  pivot_longer(cols = c(Prop_pos, Prop_neg),
               names_to = "Direction", values_to = "Proportion") %>%
  mutate(Direction = str_remove(Direction, "Prop_"),
         Model = str_replace_all(Model, x_labs))

Bigger_true_eff_tbl <- Parms_temp_bs %>%
  filter(Temp_eff == "Proportionally larger") %>%
  mutate(True_temp_eff = if_else(Correlation > 0, "Longer", "Fatter")) %>%
  janitor::tabyl(True_temp_eff)

Bigger_fatter <- Bigger_true_eff_tbl %>%
  filter(True_temp_eff == "Fatter") %>%
  mutate(percent = round(percent, 2)) %>%
  pull(percent)
Bigger_longer <- 1 - Bigger_fatter

# Inline percentages cited in the manuscript ------------------------------
pull_percent <- function(Model, Temp_eff) {
  Eval_tbl %>%
    filter(Model == {{ Model }} & Temp_eff == {{ Temp_eff }}) %>%
    mutate(per_corr = round(Prop_correct * 100, 0)) %>%
    pull(per_corr)
}

pull_percent_proportional <- function(Model, Direction, tbl = Proportional_methods_eval) {
  tbl %>%
    filter(Model == {{ Model }} & Direction == {{ Direction }}) %>%
    mutate(per_corr = round(Proportion * 100, 0)) %>%
    pull(per_corr)
}

Sli.iso_fatter_right      <- pull_percent("Sli_iso", "Fatter")
Sli.iso_longer_right      <- pull_percent("Sli_iso", "Longer")
Ryding_longer_right       <- pull_percent("Ryding",  "Longer")
Ryding_fatter_right       <- pull_percent("Ryding",  "Fatter")
Ratio_longer_right        <- pull_percent("Ratio",   "Longer")
Ratio_fatter_right        <- pull_percent("Ratio",   "Fatter")

Sli.iso_proportional_pos  <- pull_percent_proportional("SLI isometry", "pos")
Dif_sli.iso <- Sli.iso_proportional_pos - (Prop_longer * 100)
Dif_ryding  <- pull_percent_proportional("Mass as covariate", "neg") - (Prop_fatter * 100)

Sli.iso_bigger_pos <- pull_percent_proportional("SLI isometry",     "pos", tbl = Bigger_methods_eval)
Sli.est_bigger_pos <- pull_percent_proportional("SLI estimated",    "pos", tbl = Bigger_methods_eval)
Ratio_bigger_pos   <- pull_percent_proportional("Appendage / mass", "pos", tbl = Bigger_methods_eval)
Ryding_bigger_pos  <- pull_percent_proportional("Mass as covariate","pos", tbl = Bigger_methods_eval)
Bigger_longer_pct  <- round(Bigger_longer * 100, 0)

# Export --------------------------------------------------------------------
dir.create("Derived/Rds", showWarnings = FALSE)
saveRDS(
  list(
    # Settings / parameter grid (used inline in Methods prose)
    N_spp = N_spp, N_ind = N_ind,
    r_12 = r_12, r_13 = r_13, r_23 = r_23, b_avg_12 = b_avg_12,
    sma_or_ma = sma_or_ma, log_ratio = log_ratio,

    # Parameter grid + validation tables
    Parms_mat3    = Parms_mat3,
    Parms_temp_bs = Parms_temp_bs,
    Sim_fail      = Sim_fail,
    Parms_tbl3    = Parms_tbl3,
    Parms_tbl4    = Parms_tbl4,
    Mass_cor_tbl  = Mass_cor_tbl,
    x_labs        = x_labs,

    # Evaluation tables
    Eval_tbl                   = Eval_tbl,
    Proportional_methods_eval  = Proportional_methods_eval,
    Prop_fatter = Prop_fatter, Prop_longer = Prop_longer,
    Bigger_methods_eval        = Bigger_methods_eval,
    Bigger_fatter = Bigger_fatter, Bigger_longer = Bigger_longer,
    Bigger_longer_pct = Bigger_longer_pct,

    # Inline percentages
    Sli.iso_fatter_right = Sli.iso_fatter_right, Sli.iso_longer_right = Sli.iso_longer_right,
    Ryding_longer_right  = Ryding_longer_right,  Ryding_fatter_right  = Ryding_fatter_right,
    Ratio_longer_right   = Ratio_longer_right,   Ratio_fatter_right   = Ratio_fatter_right,
    Sli.iso_proportional_pos = Sli.iso_proportional_pos,
    Dif_sli.iso = Dif_sli.iso, Dif_ryding = Dif_ryding,
    Sli.iso_bigger_pos = Sli.iso_bigger_pos, Sli.est_bigger_pos = Sli.est_bigger_pos,
    Ratio_bigger_pos   = Ratio_bigger_pos,   Ryding_bigger_pos  = Ryding_bigger_pos
  ),
  "Derived/Rds/simulation_results.rds"
)
