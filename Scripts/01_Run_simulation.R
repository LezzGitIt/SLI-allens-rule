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

library(sliR)   # SLI + simulation functions; see github.com/LezzGitIt/sliR
source("Scripts/00_Key_allometry_fns.R")

# Global settings -----------------------------------------------------------
sma_or_ma <- "SMA"
log_ratio <- TRUE   # TRUE = log(A/S) and log(A²/S); FALSE = log(A)/log(S) and log(A)²/log(S)
N_ind     <- "3,000"

r_12     <- c(.3, .45, .6)
r_13     <- c(0, -.15, -.3, -.5, -.7)
r_23     <- r_13
b_avg_12 <- c(0.22, 0.33, 0.44)

# Parameter matrix ------------------------------------------------------------
# Main grid: Longer / Stouter / Proportionally smaller -- three realizations of Bergmann's rule (temperature never increases mass or appendage length).
Shape <- c("Longer", "Proportionally smaller", "Stouter")
Shape <- setNames(Shape, Shape)

Parms_mat <- expand_grid(
  b_avg_12 = b_avg_12,
  r_12 = r_12,
  r_13 = r_13,
  r_23 = r_23
) %>%
  mutate(
    # Population SMA slope of Append~Mass implied by b_avg_12 (average of OLS + SMA slope) and r_12.
    true_b_sma = 2 * b_avg_12 / (r_12 + 1),
    # Flag the "Proportionally smaller" design cells before r_13 is solved below, so the Longer/Stouter sign logic further down stays driven by the same untouched off-diagonal rows as before -- not by re-testing beta_iso_true == 0, which also happens to hold for one coincidental off-diagonal cell (b_avg_12=0.44, r_12=0.6, r_13=-0.3, r_23=-0.5; 0.3*0.55 == 0.33*0.5).
    is_proportional = r_13 == r_23,
    # Solve r_13 (instead of leaving it equal to r_23) so beta_iso_true = 0 exactly for every Proportionally smaller cell regardless of allometric scaling category -- r_13 = r_23 only zeroed it for Hypoallometric species; see Project_notes.md / Ground_truth_explainer.qmd "Round 7" for the derivation. Off-diagonal (Longer/Stouter) rows keep their original r_13.
    r_13 = if_else(is_proportional, 0.33 * r_23 / true_b_sma, r_13),
    # Closed-form sign of SLI-isometry's population coefficient on temperature -- this simulation's ground truth for "true" shape change (see Project_notes.md for derivation).
    beta_iso_true = r_13 * true_b_sma - 0.33 * r_23,
    Temp_eff = case_when(
      is_proportional ~ Shape[2],
      beta_iso_true > 0 ~ Shape[1],
      beta_iso_true < 0 ~ Shape[3]
    ),
    Strength = abs(r_13 - r_23),
    is_proportional = NULL
  ) %>%
  mutate(Scaling = case_when(
    b_avg_12 < 0.33 ~ "Hypoallometry",
    b_avg_12 > 0.33 ~ "Hyperallometry",
    near(b_avg_12, 0.33) ~ "Isometry"
  ))

# Proportionally larger: mirror of the negative "Proportionally smaller" diagonal, on the positive-r_23 side -- temperature increases wing and mass with equal, not differential, pull. r_13 is solved (not set equal to r_23) so beta_iso_true = 0 exactly, same fix as Proportionally smaller above. Simulated analogue of empirical "Inverse Bergmann's" species.
r_temp_pos <- c(.15, .3, .5, .7)

Parms_big <- expand_grid(
  b_avg_12 = b_avg_12,
  r_12     = r_12,
  r_23     = r_temp_pos
) %>%
  mutate(true_b_sma    = 2 * b_avg_12 / (r_12 + 1),
         r_13          = 0.33 * r_23 / true_b_sma,
         beta_iso_true = r_13 * true_b_sma - 0.33 * r_23,
         Temp_eff      = "Proportionally larger",
         Strength      = abs(r_13 - r_23),
         Scaling       = case_when(
           b_avg_12 < 0.33 ~ "Hypoallometry",
           b_avg_12 > 0.33 ~ "Hyperallometry",
           near(b_avg_12, 0.33) ~ "Isometry"
         ))

Parms_mat2 <- bind_rows(Parms_mat, Parms_big) %>%
  mutate(
    Temp_eff = factor(Temp_eff,
                      levels = c("Stouter", "Proportionally smaller", "Proportionally larger", "Longer")),
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
    Temp_eff == "Stouter"  & cor_allometry_int > -.2 |
    Temp_eff == "Stouter"  & cor_allometry     > -.2 |
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
           # FALSE branch fixed: this must be log(A^2)/log(S) = 2*log(A)/log(S), the ratio-of-logs
           # analog of squaring A before dividing -- Append_log^2/Mass_log (squaring the already-
           # logged value) was a different, inconsistent quantity with no such interpretation.
           Append2_mass = if (log_ratio) 2*Append_log - Mass_log  else 2*Append_log / Mass_log,
           # Ratio-of-logs robustness check (supplement only -- see Parms_tbl4_altratio below):
           # log(A)/log(S) and log(A^2)/log(S), computed unconditionally alongside the primary
           # log_ratio-gated columns above so both constructions are always available from one run.
           Append_mass_alt  = Append_log / Mass_log,
           Append2_mass_alt = 2 * Append_log / Mass_log) %>%
    sliR::calc_sli(b_sli = 0.33,    rename_col = "sli_isometry") %>%
    sliR::calc_sli(b_sli = est_b_sma, rename_col = "sli_estimated")
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
  # Ratio-of-logs robustness check (supplement only): see Append_mass_alt/Append2_mass_alt above.
  ratio_alt_app  <- lm(Append_mass_alt  ~ Temp_inc, data = Sim_df_s, na.action = na.exclude)
  ratio2_alt_app <- lm(Append2_mass_alt ~ Temp_inc, data = Sim_df_s, na.action = na.exclude)

  # 95% CI on the Temp_inc coefficient for each method -- used downstream to define
  # a stricter correctness criterion (point estimate right-signed AND CI excludes
  # zero), not just point-estimate sign. "_app" suffix kept on every one of the six
  # (Ratio/Ratio2 included) so a single pivot_longer(.value) can reshape coef/ci_lo/ci_hi
  # together below.
  ci_temp_inc <- function(mod) confint(mod)["Temp_inc", ]

  tibble(
    coef_sli_iso_app   = coef(sli_iso_app)["Temp_inc"],
    coef_sli_est_app   = coef(sli_est_app)["Temp_inc"],
    coef_ols_resid_app = coef(ols_resid_app)["Temp_inc"],
    coef_ryding_app    = coef(ryding_app)["Temp_inc"],
    coef_ratio_app     = coef(ratio_app)["Temp_inc"],
    coef_ratio2_app    = coef(ratio2_app)["Temp_inc"],
    ci_lo_sli_iso_app   = ci_temp_inc(sli_iso_app)[1],
    ci_hi_sli_iso_app   = ci_temp_inc(sli_iso_app)[2],
    ci_lo_sli_est_app   = ci_temp_inc(sli_est_app)[1],
    ci_hi_sli_est_app   = ci_temp_inc(sli_est_app)[2],
    ci_lo_ols_resid_app = ci_temp_inc(ols_resid_app)[1],
    ci_hi_ols_resid_app = ci_temp_inc(ols_resid_app)[2],
    ci_lo_ryding_app    = ci_temp_inc(ryding_app)[1],
    ci_hi_ryding_app    = ci_temp_inc(ryding_app)[2],
    ci_lo_ratio_app     = ci_temp_inc(ratio_app)[1],
    ci_hi_ratio_app     = ci_temp_inc(ratio_app)[2],
    ci_lo_ratio2_app    = ci_temp_inc(ratio2_app)[1],
    ci_hi_ratio2_app    = ci_temp_inc(ratio2_app)[2],
    # Ratio-of-logs robustness check (supplement only): coef/CI columns follow the same
    # "_app" naming convention on purpose, so Parms_tbl4's pivot_longer() below sweeps
    # them in as extra Model levels automatically -- split back out immediately after
    # (Parms_tbl4_altratio) so the exported Parms_tbl4 that feeds main-text figures is
    # unaffected.
    coef_ratio_alt_app  = coef(ratio_alt_app)["Temp_inc"],
    coef_ratio2_alt_app = coef(ratio2_alt_app)["Temp_inc"],
    ci_lo_ratio_alt_app  = ci_temp_inc(ratio_alt_app)[1],
    ci_hi_ratio_alt_app  = ci_temp_inc(ratio_alt_app)[2],
    ci_lo_ratio2_alt_app = ci_temp_inc(ratio2_alt_app)[1],
    ci_hi_ratio2_alt_app = ci_temp_inc(ratio2_alt_app)[2],
    est_b_sma          = coefs$est_b_sma,
    est_b_ols          = coefs$est_b_ols,
    # Per-species Pearson correlation between each method's individual-level
    # metric and body mass (@fig-mass-cor-sim in supplementary_info.qmd),
    # confirming in simulated data the mechanical ratio/mass dependence shown
    # empirically in the main text (@tbl-ratio-mass-summary). Multiple regression
    # is excluded: it has no individual-level shape metric to correlate, only
    # a model coefficient.
    cor_ratio          = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$Append_mass)),
    cor_ratio2         = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$Append2_mass)),
    cor_ols_resid      = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$resid_ols)),
    cor_sli_iso        = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$sli_isometry)),
    cor_sli_est        = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$sli_estimated)),
    # Ratio-of-logs robustness check (supplement only)
    cor_ratio_alt      = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$Append_mass_alt)),
    cor_ratio2_alt     = as.numeric(cor(Sim_df_s$Mass_log, Sim_df_s$Append2_mass_alt))
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
  "Ryding"    = "Multiple regression",
  "Ratio2"    = "Appendage² / mass",
  "Ratio"     = "Appendage / mass",
  # Ratio-of-logs robustness check (supplement only)
  "Ratio_alt"  = "Appendage / mass (log-of-logs)",
  "Ratio2_alt" = "Appendage² / mass (log-of-logs)"
)

# Est_correct requires the whole 95% CI, not just the point estimate, to fall on the correct
# side of zero -- a species simulated as Stouter with a negative point estimate but a CI spanning
# zero does not count as correctly classified. This necessarily implies the point-estimate sign
# is also correct (the point estimate always lies within its own CI). For Proportionally
# smaller/larger species, the true effect is exactly 0 by construction (Effect classification),
# so the criterion reverses: a method is correct there if its CI *includes* zero (correctly
# failing to detect a nonexistent effect), not if it excludes it. Computed here (not per-qmd) so
# every downstream consumer of Parms_tbl4 -- Eval_tbl, both fig-eval-style calibration scatters,
# Results-paragraph percentages -- shares one definition of "correct."
Parms_tbl4_full <- Parms_tbl3 %>%
  # Drop the alt-ratio mass-correlation columns before pivoting -- they aren't matched by
  # the coef/ci_lo/ci_hi "_app" pattern below, so left in they'd otherwise ride along as
  # extra id_cols on every row of Parms_tbl4, including the original six methods'. Keeping
  # Parms_tbl4's schema exactly as before is what lets the byte-identical check above pass
  # cleanly, not just "functionally equivalent with harmless extra columns."
  dplyr::select(-cor_ratio_alt, -cor_ratio2_alt) %>%
  pivot_longer(
    cols = matches("^(coef|ci_lo|ci_hi)_.+_app$"),
    names_to = c(".value", "Model"),
    names_pattern = "(coef|ci_lo|ci_hi)_(.+)_app"
  ) %>%
  rename(b_temp_inc = coef) %>%
  mutate(Model = str_to_sentence(Model)) %>%
  mutate(b_dir = if_else(b_temp_inc < 0, "Neg", "Pos")) %>%
  mutate(Est_correct = case_when(
    Temp_eff == "Longer" ~ ci_lo > 0,
    Temp_eff == "Stouter" ~ ci_hi < 0,
    Temp_eff %in% c("Proportionally smaller", "Proportionally larger") ~ ci_lo <= 0 & ci_hi >= 0,
    .default = FALSE
  ), .by = Model)

# Ratio-of-logs robustness check (supplement only): Ratio_alt/Ratio2_alt (log(A)/log(S) and
# log(A^2)/log(S)) got swept into Parms_tbl4_full above by the same generic pivot_longer() as
# the primary six methods, since their coef/CI columns follow the same "_app" naming convention
# (see extract_coefs() above). Split them back out immediately, BEFORE any downstream object is
# built, so every main-text consumer of "Parms_tbl4" -- fig-compare-approaches and fig-eval in
# Allens_methods_sim.qmd (both read Parms_tbl4 with no Model filter), and
# supplementary_info.qmd's Parms_tbl4_bigger -- sees exactly the original six methods, unchanged.
main_models <- c("Sli_iso", "Sli_est", "Ols_resid", "Ryding", "Ratio", "Ratio2")
Parms_tbl4_altratio <- Parms_tbl4_full %>% filter(Model %in% c(main_models, "Ratio_alt", "Ratio2_alt"))
Parms_tbl4          <- Parms_tbl4_full %>% filter(Model %in% main_models)

# Per-species correlation between each method's metric and body mass ---------
Mass_cor_tbl <- Parms_tbl3 %>%
  # Drop the ratio-of-logs robustness-check columns first -- same reason as
  # Parms_tbl4_full above: pivot_longer()'s unlisted columns become id_cols by default,
  # and Mass_cor_tbl_altratio (below) already covers this comparison from Parms_tbl3
  # directly, so nothing is lost by excluding them here.
  dplyr::select(-matches("_alt")) %>%
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
Eval_tbl <- Parms_tbl4 %>%
  summarize(Prop_correct = sum(Est_correct, na.rm = TRUE) / n(),
            .by = c(Model, Temp_eff)) %>%
  mutate(Prop_correct = round(Prop_correct, 2)) %>%
  arrange(Model)

# Ratio-of-logs robustness check (supplement only): does the log(A)/log(S) /
# log(A^2)/log(S) construction change the overall accuracy trends versus the primary
# log(A/S) / log(A^2/S) ratios, or their mass-dependence? Mirrors Eval_tbl/Mass_cor_tbl
# above exactly, built from the Parms_tbl4_altratio split-off (Ratio_alt/Ratio2_alt
# plus the four non-ratio methods for reference) rather than the main Parms_tbl4.
Eval_tbl_altratio <- Parms_tbl4_altratio %>%
  summarize(Prop_correct = sum(Est_correct, na.rm = TRUE) / n(),
            .by = c(Model, Temp_eff)) %>%
  mutate(Prop_correct = round(Prop_correct, 2)) %>%
  arrange(Model)

Mass_cor_tbl_altratio <- Parms_tbl3 %>%
  pivot_longer(
    cols = c(cor_ratio, cor_ratio2, cor_ratio_alt, cor_ratio2_alt, cor_ols_resid, cor_sli_iso, cor_sli_est),
    names_to  = "Model",
    values_to = "r_mass"
  ) %>%
  mutate(Model = str_remove_all(Model, "cor_"),
         Model = str_to_sentence(Model))

# OLS-anchored robustness check (supplementary): recompute the ground truth with 0.33 compared against an OLS-scale reference slope instead of the adopted SMA-scale one, and re-score all six methods' sign-agreement against it. See Project_notes.md / Ground_truth_explainer.qmd "Round 3" for the derivation and the SMA-vs-OLS rationale it is checking (main text, Approach classification). Restricted to Longer/Stouter, matching Eval_tbl's convention: beta_iso_true is exactly 0 for both Proportional categories under the SMA anchor by design (Effect classification), so a signed sign-agreement comparison isn't meaningful for them under that anchor.
Ols_anchor_df <- Parms_tbl4 %>%
  filter(Temp_eff %in% c("Longer", "Stouter")) %>%
  mutate(Model = str_replace_all(Model, x_labs),
         b_true_ols        = r_12 * true_b_sma,
         beta_iso_true_ols = r_13 * b_true_ols - 0.33 * r_23)

# "Correct" requires the whole 95% CI, not just the point estimate, on the correct side of
# zero -- same criterion as Est_correct above, applied here to both candidate ground truths.
ci_correct <- function(true_val, ci_lo, ci_hi) (true_val > 0 & ci_lo > 0) | (true_val < 0 & ci_hi < 0)

Ols_anchor_summary_tbl <- Ols_anchor_df %>%
  group_by(Model) %>%
  summarize(
    `r (SMA-anchored)` = round(cor(beta_iso_true, b_temp_inc), 2),
    `% correct (SMA)`  = round(mean(ci_correct(beta_iso_true, ci_lo, ci_hi)) * 100, 1),
    `r (OLS-anchored)` = round(cor(beta_iso_true_ols, b_temp_inc), 2),
    `% correct (OLS)`  = round(mean(ci_correct(beta_iso_true_ols, ci_lo, ci_hi)) * 100, 1),
    .groups = "drop"
  )

# Figure-only companion to Ols_anchor_df: widens the OLS-anchored comparison to also include
# Proportionally smaller species, matching main text's fig-eval/Eval_scatter_df convention
# (excludes only Proportionally larger, the out-of-main-grid scenario). This can't just reuse
# Ols_anchor_df/Ols_anchor_summary_tbl above, which stay Longer/Stouter-only so Table 1's
# figures in the surrounding prose are unaffected by widening the figure's scope. Unlike the SMA
# anchor, which pins Proportionally smaller species to an exactly-null true effect by construction
# (main text, Effect classification) and so needs Est_correct's reversed "CI includes zero" rule,
# beta_iso_true_ols is not pinned to zero for these species -- the SMA-vs-OLS anchor choice is
# exactly what breaks that exact nullity -- so the same sign-based ci_correct() used for
# Longer/Stouter applies here without a special case.
Ols_anchor_fig_df <- Parms_tbl4 %>%
  filter(Temp_eff != "Proportionally larger") %>%
  mutate(Model = str_replace_all(Model, x_labs),
         b_true_ols        = r_12 * true_b_sma,
         beta_iso_true_ols = r_13 * b_true_ols - 0.33 * r_23,
         Est_correct_ols   = ci_correct(beta_iso_true_ols, ci_lo, ci_hi))

# Precomputed r/pct_correct pair for the fig-eval-style calibration scatter in supplementary_info.qmd -- computed here, not in the qmd, matching this section's "no new stats inside the qmd" convention.
Ols_anchor_stats <- Ols_anchor_fig_df %>%
  group_by(Model) %>%
  summarize(r = round(cor(beta_iso_true_ols, b_temp_inc), 2),
            pct_correct = round(mean(Est_correct_ols) * 100, 1),
            .groups = "drop")

pull_ols_pct <- function(model, col, tbl = Ols_anchor_summary_tbl) {
  tbl %>% filter(Model == {{ model }}) %>% pull({{ col }})
}

Sliiso_pct_sma   <- pull_ols_pct("SLI isometry",      `% correct (SMA)`)
Sliiso_pct_ols   <- pull_ols_pct("SLI isometry",      `% correct (OLS)`)
Ratio_pct_sma    <- pull_ols_pct("Appendage / mass",  `% correct (SMA)`)
Ratio_pct_ols    <- pull_ols_pct("Appendage / mass",  `% correct (OLS)`)
Ryding_pct_sma   <- pull_ols_pct("Multiple regression", `% correct (SMA)`)
Ryding_pct_ols   <- pull_ols_pct("Multiple regression", `% correct (OLS)`)
Olsresid_pct_sma <- pull_ols_pct("OLS residuals",     `% correct (SMA)`)
Olsresid_pct_ols <- pull_ols_pct("OLS residuals",     `% correct (OLS)`)

# SLI-isometry / beta_T rescaling check (supplementary): beta_iso_true only captures
# Cov(log SLI, Temp), while SLI-isometry's actual fitted coefficient is that covariance divided by
# sd(log SLI) too (a correlation, not a covariance) -- a species-specific variance term beta_iso_true
# never accounted for. V = Var(Append - 0.33*Mass)/sd_M^2 = true_b_sma^2 + 0.33^2 - 2*0.33*r_12*true_b_sma;
# dividing beta_iso_true by sqrt(V) restores the missing term and should put SLI-isometry almost
# exactly on the 1:1 line (r ~ 0.9999) against a raw r ~ 0.97 -- see Project_notes.md /
# Ground_truth_explainer.qmd for the derivation. Scoped to Longer/Stouter/Proportionally smaller,
# matching fig-eval; Proportionally smaller collapses to 0 under both the raw and rescaled versions.
Sli_iso_rescale_df <- Parms_tbl4 %>%
  filter(Temp_eff != "Proportionally larger", Model == "Sli_iso") %>%
  mutate(V = true_b_sma^2 + 0.33^2 - 2 * 0.33 * r_12 * true_b_sma,
         beta_T_rescaled = beta_iso_true / sqrt(V))

Sli_iso_rescale_stats <- Sli_iso_rescale_df %>%
  summarize(r_raw = round(cor(beta_iso_true, b_temp_inc), 4),
            r_rescaled = round(cor(beta_T_rescaled, b_temp_inc), 4))
r_sli_iso_raw      <- Sli_iso_rescale_stats$r_raw
r_sli_iso_rescaled <- Sli_iso_rescale_stats$r_rescaled

Proportional_methods_eval <- Parms_tbl4 %>%
  filter(Temp_eff == "Proportionally smaller") %>%
  summarize(Prop_pos = sum(b_dir == "Pos") / n(), .by = c(Model, Temp_eff)) %>%
  mutate(Prop_neg = 1 - Prop_pos) %>%
  pivot_longer(cols = c(Prop_pos, Prop_neg),
               names_to = "Direction", values_to = "Proportion") %>%
  mutate(Direction = str_remove(Direction, "Prop_"),
         Model = str_replace_all(Model, x_labs))

# True benchmark is now an exact analytic constant, not a quantity tabulated from beta_iso_true's sign (which is 0, not >0 or <0, for every row in both categories now that r_13 is solved above) -- an unbiased method's point-estimate sign should split ~50/50 by sampling noise alone, since the true effect is exactly null in every species in both categories.
stopifnot(max(abs(Parms_mat3$beta_iso_true[Parms_mat3$Temp_eff == "Proportionally smaller"])) < 1e-9)
stopifnot(max(abs(Parms_mat3$beta_iso_true[Parms_mat3$Temp_eff == "Proportionally larger"])) < 1e-9)

Prop_true_null_pct <- 50
Prop_fatter   <- 0.5
Prop_longer   <- 0.5

# "Proportionally larger" scenario: same "no true shape signal" evaluation as Proportionally smaller, just mirrored to the positive-r_13/r_23 diagonal (species getting bigger, not smaller). Kept as a parallel/duplicate computation rather than generalizing Proportional_methods_eval, so those objects remain unaffected.
Bigger_methods_eval <- Parms_tbl4 %>%
  filter(Temp_eff == "Proportionally larger") %>%
  summarize(Prop_pos = sum(b_dir == "Pos") / n(), .by = c(Model, Temp_eff)) %>%
  mutate(Prop_neg = 1 - Prop_pos) %>%
  pivot_longer(cols = c(Prop_pos, Prop_neg),
               names_to = "Direction", values_to = "Proportion") %>%
  mutate(Direction = str_remove(Direction, "Prop_"),
         Model = str_replace_all(Model, x_labs))

Bigger_fatter <- 0.5
Bigger_longer <- 0.5

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

Sli.iso_fatter_right      <- pull_percent("Sli_iso", "Stouter")
Sli.iso_longer_right      <- pull_percent("Sli_iso", "Longer")
Ryding_longer_right       <- pull_percent("Ryding",  "Longer")
Ryding_fatter_right       <- pull_percent("Ryding",  "Stouter")
Ratio_longer_right        <- pull_percent("Ratio",   "Longer")
Ratio_fatter_right        <- pull_percent("Ratio",   "Stouter")
Sli.est_longer_right      <- pull_percent("Sli_est", "Longer")
Sli.est_fatter_right      <- pull_percent("Sli_est", "Stouter")

Ratio2_fatter_right       <- pull_percent("Ratio2", "Stouter")
# "Wrong" here matches Est_correct's CI-based criterion (Evaluation section above), not just
# a wrong-signed point estimate -- a right-signed point estimate whose CI includes zero also counts as wrong.
Ratio2_fatter_wrong_n <- Parms_tbl4 %>%
  filter(Temp_eff == "Stouter", Model == "Ratio2", !Est_correct) %>% nrow()
Fatter_n <- Parms_tbl4 %>% filter(Temp_eff == "Stouter", Model == "Ratio2") %>% nrow()

Sli.iso_proportional_pos  <- pull_percent_proportional("SLI isometry", "pos")
Sli.est_proportional_pos  <- pull_percent_proportional("SLI estimated", "pos")
Ratio_proportional_pos    <- pull_percent_proportional("Appendage / mass", "pos")
Ratio2_proportional_fatter <- pull_percent_proportional("Appendage² / mass", "neg")
Dif_sli.iso <- Sli.iso_proportional_pos - Prop_true_null_pct
Dif_sli.est <- Sli.est_proportional_pos - Prop_true_null_pct
Dif_ryding  <- pull_percent_proportional("Multiple regression", "neg") - Prop_true_null_pct
Dif_ratio   <- Ratio_proportional_pos - Prop_true_null_pct

Sli.iso_bigger_pos <- pull_percent_proportional("SLI isometry",     "pos", tbl = Bigger_methods_eval)
Sli.est_bigger_pos <- pull_percent_proportional("SLI estimated",    "pos", tbl = Bigger_methods_eval)
Ratio_bigger_pos   <- pull_percent_proportional("Appendage / mass", "pos", tbl = Bigger_methods_eval)
Ryding_bigger_pos  <- pull_percent_proportional("Multiple regression","pos", tbl = Bigger_methods_eval)
Bigger_longer_pct  <- round(Bigger_longer * 100, 0)

# CI-based null-detection accuracy: % of species where the 95% CI correctly *includes* zero
# (Est_correct's reversed criterion for the two Proportional categories, Evaluation section above).
# Eval_tbl already carries these rows for every Model, so pull_percent() works unmodified.
Sli.iso_proportional_null_correct   <- pull_percent("Sli_iso",   "Proportionally smaller")
Sli.est_proportional_null_correct   <- pull_percent("Sli_est",   "Proportionally smaller")
Ratio_proportional_null_correct     <- pull_percent("Ratio",     "Proportionally smaller")
Ratio2_proportional_null_correct    <- pull_percent("Ratio2",    "Proportionally smaller")
Ryding_proportional_null_correct    <- pull_percent("Ryding",    "Proportionally smaller")
Olsresid_proportional_null_correct  <- pull_percent("Ols_resid", "Proportionally smaller")

Sli.iso_bigger_null_correct   <- pull_percent("Sli_iso",   "Proportionally larger")
Sli.est_bigger_null_correct   <- pull_percent("Sli_est",   "Proportionally larger")
Ratio_bigger_null_correct     <- pull_percent("Ratio",     "Proportionally larger")
Ratio2_bigger_null_correct    <- pull_percent("Ratio2",    "Proportionally larger")
Ryding_bigger_null_correct    <- pull_percent("Ryding",    "Proportionally larger")
Olsresid_bigger_null_correct  <- pull_percent("Ols_resid", "Proportionally larger")

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

    # Ratio-of-logs robustness check (supplement only)
    Parms_tbl4_altratio   = Parms_tbl4_altratio,
    Eval_tbl_altratio      = Eval_tbl_altratio,
    Mass_cor_tbl_altratio  = Mass_cor_tbl_altratio,

    # Evaluation tables
    Eval_tbl                   = Eval_tbl,
    Proportional_methods_eval  = Proportional_methods_eval,
    Prop_fatter = Prop_fatter, Prop_longer = Prop_longer,
    Prop_true_null_pct = Prop_true_null_pct,
    Bigger_methods_eval        = Bigger_methods_eval,
    Bigger_fatter = Bigger_fatter, Bigger_longer = Bigger_longer,
    Bigger_longer_pct = Bigger_longer_pct,

    # OLS-anchored robustness check (supplementary)
    Ols_anchor_df          = Ols_anchor_df,
    Ols_anchor_fig_df      = Ols_anchor_fig_df,
    Ols_anchor_summary_tbl = Ols_anchor_summary_tbl,
    Ols_anchor_stats       = Ols_anchor_stats,
    Sliiso_pct_sma = Sliiso_pct_sma, Sliiso_pct_ols = Sliiso_pct_ols,
    Ratio_pct_sma  = Ratio_pct_sma,  Ratio_pct_ols  = Ratio_pct_ols,
    Ryding_pct_sma = Ryding_pct_sma, Ryding_pct_ols = Ryding_pct_ols,
    Olsresid_pct_sma = Olsresid_pct_sma, Olsresid_pct_ols = Olsresid_pct_ols,

    # Inline percentages
    Sli.iso_fatter_right = Sli.iso_fatter_right, Sli.iso_longer_right = Sli.iso_longer_right,
    Ryding_longer_right  = Ryding_longer_right,  Ryding_fatter_right  = Ryding_fatter_right,
    Ratio_longer_right   = Ratio_longer_right,   Ratio_fatter_right   = Ratio_fatter_right,
    Sli.est_longer_right = Sli.est_longer_right, Sli.est_fatter_right = Sli.est_fatter_right,
    Ratio2_fatter_right  = Ratio2_fatter_right,  Ratio2_fatter_wrong_n = Ratio2_fatter_wrong_n,
    Fatter_n = Fatter_n,
    Sli.iso_proportional_pos = Sli.iso_proportional_pos,
    Sli.est_proportional_pos = Sli.est_proportional_pos,
    Ratio_proportional_pos = Ratio_proportional_pos,
    Ratio2_proportional_fatter = Ratio2_proportional_fatter,
    Dif_sli.iso = Dif_sli.iso, Dif_sli.est = Dif_sli.est, Dif_ryding = Dif_ryding,
    Dif_ratio = Dif_ratio,
    Sli.iso_bigger_pos = Sli.iso_bigger_pos, Sli.est_bigger_pos = Sli.est_bigger_pos,
    Ratio_bigger_pos   = Ratio_bigger_pos,   Ryding_bigger_pos  = Ryding_bigger_pos,

    # CI-based null-detection accuracy (Proportionally smaller/larger)
    Sli.iso_proportional_null_correct  = Sli.iso_proportional_null_correct,
    Sli.est_proportional_null_correct  = Sli.est_proportional_null_correct,
    Ratio_proportional_null_correct    = Ratio_proportional_null_correct,
    Ratio2_proportional_null_correct   = Ratio2_proportional_null_correct,
    Ryding_proportional_null_correct   = Ryding_proportional_null_correct,
    Olsresid_proportional_null_correct = Olsresid_proportional_null_correct,
    Sli.iso_bigger_null_correct  = Sli.iso_bigger_null_correct,
    Sli.est_bigger_null_correct  = Sli.est_bigger_null_correct,
    Ratio_bigger_null_correct    = Ratio_bigger_null_correct,
    Ratio2_bigger_null_correct   = Ratio2_bigger_null_correct,
    Ryding_bigger_null_correct   = Ryding_bigger_null_correct,
    Olsresid_bigger_null_correct = Olsresid_bigger_null_correct,

    # SLI-isometry variance-rescaling check (supplementary)
    Sli_iso_rescale_df = Sli_iso_rescale_df,
    r_sli_iso_raw      = r_sli_iso_raw,
    r_sli_iso_rescaled = r_sli_iso_rescaled
  ),
  "Derived/Rds/simulation_results.rds"
)
