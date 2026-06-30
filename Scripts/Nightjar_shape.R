## Analysis -- Shapeshifting in nightjars ##


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(smatr)
library(cowplot)
library(broom)
library(ggpmisc)
library(lavaan)
ggplot2::theme_set(theme_cowplot())

source("Scripts/Key_allometry_fns.R")
source("Scripts/Key_causal_fns.R")
nj_raw <- read.csv("Data/Capri_BA_compare03.29.26.csv")

# Run parameters -----------------------------------------------------------
include_sem     <- TRUE   # include lavaan SEM in the analysis
                           # TRUE:  all methods use Wing+Mass+Tail complete cases (comparability)
                           # FALSE: use Wing+Mass complete cases for OLS/SMA/SLI/ratio only
control_age_sex <- TRUE   # include Age/Sex as covariates in the SEM
                           # Sex: all species; Age: Nightjar and Whip-poor-will only

# Quick data summary -------------------------------------------------------
nj_raw %>%
  mutate(log_wing = log(Wing.comb), log_mass = log(Mass.comb)) %>%
  summarize(
    mn.wing = mean(log_wing, na.rm = TRUE),
    sd.wing = sd(log_wing, na.rm = TRUE),
    mn.mass = mean(log_mass, na.rm = TRUE),
    sd.mass = sd(log_mass, na.rm = TRUE),
    b_sma   = sd.wing / sd.mass,
    .by = Species)

# Formatting --------------------------------------------------------------
# NOTE on sign convention: B.Lat (latitude, degrees N) is the temperature
# proxy — higher latitude is colder. A NEGATIVE B.Lat coefficient therefore
# indicates Allen's rule (longer appendages at lower/warmer latitudes).
# This is the reverse of the simulation convention where higher Temp = warmer.

nj_df <- nj_raw %>%
  dplyr::select(Species,
                Wing = Wing.comb, Mass = Mass.comb, Tail = Tail.comb,
                B.Lat, Age, Sex) %>%
  drop_na(Wing, Mass) %>%
  mutate(log_wing   = log(Wing),
         log_mass   = log(Mass),
         log_tail   = log(Tail),
         wing_mass  = Wing / Mass,
         wing2_mass = Wing^2 / Mass)

# Visualize regression approaches -----------------------------------------
# Compare three line-fitting methods
nj_df %>%
  ggplot(aes(x = log_mass, y = log_wing)) +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, color = "red") +
  ggpmisc::stat_ma_line(method = "MA", se = FALSE, color = "orange") +
  ggpmisc::stat_ma_line(method = "SMA", linetype = "dotted", se = FALSE, color = "blue") +
  facet_wrap(~Species)

# Examine variance in X & Y
nj_df %>% group_by(Species) %>%
  summarise(var_mass = var(log_mass), var_wing = var(log_wing)) %>%
  mutate(ratio = sqrt(var_wing / var_mass))

# Create lists by Species -------------------------------------------------
# List 1: Wing + Mass complete (full dataset; baseline when include_sem = FALSE)
nj_df_l <- nj_df %>% group_split(Species)
names(nj_df_l) <- c("Nighthawk", "Nightjar", "Whip-poor-will")

# List 2: Wing + Mass + Tail complete (required for 3-indicator SEM)
nj_df_sem_l <- nj_df %>%
  drop_na(Tail) %>%
  group_split(Species)
names(nj_df_sem_l) <- c("Nighthawk", "Nightjar", "Whip-poor-will")

# Analysis dataset: when include_sem = TRUE all methods use the Tail-complete
# dataset so estimates are directly comparable across approaches.
nj_df_analysis_l <- if (include_sem) nj_df_sem_l else nj_df_l

# Allometric correlation --------------------------------------------------
# SMA is appropriate only when there is a meaningful mass-wing correlation
# (Smith 2009). Sample sizes reflect post-morphology-filter N for the chosen
# analysis dataset.
Spp_metadata <- map(nj_df_analysis_l, \(df){
  b_sma <- sd(df$log_wing, na.rm = TRUE) / sd(df$log_mass, na.rm = TRUE)
  tibble(num_obs = nrow(df),
         b_sma   = b_sma,
         b_ols   = cor(df$log_wing, df$log_mass, use = "complete.obs") * b_sma)
}) %>% list_rbind(names_to = "Species")

Cors_tbl <- map(nj_df_analysis_l, \(df){
  tidy(lm(Mass ~ Wing, data = df))
}) %>% list_rbind(names_to = "Species") %>%
  filter(term == "Wing") %>%
  dplyr::select(Species, estimate, p.value) %>%
  rename(cor_mw = estimate, p_mw = p.value)

Spp_metadata <- left_join(Spp_metadata, Cors_tbl, by = "Species")
Spp_metadata

# Coefficient of variation ------------------------------------------------
# Very low CVs, suggesting that the variance in wing is very low compared to the variance in mass (relative to the means)
nj_df %>% summarize(N = n(),
                    lambda = calc_lambda(x = Mass, y = Wing),
                    .by = Species) %>%
  mutate(lambda = round(lambda, 2))
0.0296^2 / .102^2 # Example, European nightjar

# Test assumptions --------------------------------------------------------
# OLS assumptions
ols_mod_l <- map(nj_df_analysis_l, \(df){
  lm(log_wing ~ log_mass + B.Lat, data = df)
})

# Some departure from homoskedasticity
map(ols_mod_l, \(ols_mod){
  plot(ols_mod, which = 1)
  plot(ols_mod, which = 2)
})

# SMA assumptions
sma_mod_l <- map(nj_df_analysis_l, \(df){
  sma(log_wing ~ log_mass, data = df, method = "SMA")
})

map(sma_mod_l, \(sma_mod){
  plot(sma_mod, which = "residual")
  plot(sma_mod, which = "qq")
})

# Prep data --------------------------------------------------------------
## NOTE: if you scale first, then the variance of both log_mass & log_wing is 1,
## & SMA slope = MA slope = 1 × OLS slope so these are identical.

nj_df_l2 <- map(nj_df_analysis_l, \(df){
  ols_mod   <- lm(log_wing ~ log_mass, data = df)
  sma_mod   <- sma(log_wing ~ log_mass, data = df, method = "SMA")
  est_b_sma <- coef(sma_mod)['slope']
  df %>%
    mutate(resid_ols = residuals(ols_mod),
           resid_sma = residuals(sma_mod)) %>%
    calc_sli(b_sli = 0.33,      Append = Wing, rename_col = "sli_isometry") %>%
    calc_sli(b_sli = est_b_sma, Append = Wing, rename_col = "sli_estimated")
})

# Scale by species
nj_df_l3 <- map(nj_df_l2, \(df){
  df %>% mutate(across(where(is.numeric), scale))
})

# Inspect correlations with body size (mass in this case)
map(nj_df_l3, \(df){
  df %>% summarize(wm_m   = cor(wing_mass, Mass),
                   w2m_m  = cor(wing2_mass, Mass),
                   resid_m = cor(resid_sma, Mass))
})


# Run models & extract parms ----------------------------------------------
parms_df <- map(nj_df_l3, \(df){
  mod_resid_ols   <- lm(resid_ols   ~ B.Lat,        data = df) %>% tidy() %>% mutate(Approach = "Resid_ols")
  mod_coef_ols    <- lm(Wing        ~ Mass + B.Lat,  data = df) %>% tidy() %>% mutate(Approach = "Ryding")
  mod_coef_ratio  <- lm(wing_mass   ~ B.Lat,         data = df) %>% tidy() %>% mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass  ~ B.Lat,         data = df) %>% tidy() %>% mutate(Approach = "Ratio2")
  mod_sli_iso     <- lm(sli_isometry  ~ B.Lat,       data = df) %>% tidy() %>% mutate(Approach = "Sli_iso")
  mod_sli_est     <- lm(sli_estimated ~ B.Lat,       data = df) %>% tidy() %>% mutate(Approach = "Sli_est")
  bind_rows(mod_sli_iso, mod_sli_est, mod_resid_ols, mod_coef_ols, mod_coef_ratio, mod_coef_ratio2)
}) %>% list_rbind(names_to = "Species") %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

# SEM with 3 indicators: Mass (anchor) + Wing + Tail ---------------------
# Only computed when include_sem = TRUE.
# fit_lavaan_sem() requires pre-z-scored data; data are z-scored per species
# before calling. Age/Sex filters applied where control_age_sex = TRUE.

if (include_sem) {

  # Per-species covariate helper (Sex: all; Age: Nightjar & Whip-poor-will)
  get_sem_covs <- function(sp) {
    if (!control_age_sex) return(character(0))
    covs <- "Sex"
    if (sp %in% c("Nightjar", "Whip-poor-will")) covs <- c(covs, "Age")
    covs
  }

  nj_df_sem_l2 <- imap(nj_df_analysis_l, \(df, sp){
    out <- df %>% mutate(across(where(is.numeric), scale))
    covs <- get_sem_covs(sp)
    if ("Sex" %in% covs) out <- out %>% filter(!is.na(Sex) & Sex != "U")
    if ("Age" %in% covs) out <- out %>% filter(!is.na(Age) & Age != "Unknown")
    out
  })

  sem_parms <- imap(nj_df_sem_l2, \(df, sp) {
    res <- fit_lavaan_sem(df,
      mass_name    = "log_mass",
      append_names = c("log_wing", "log_tail"),
      temp_name    = "B.Lat",
      labels       = c("wing", "tail"),
      size_covs    = get_sem_covs(sp))
    tibble(
      term      = "B.Lat",
      estimate  = c(res$coef_sem_wing,      res$coef_sem_tail),
      std.error = c(res$se_sem_wing,        res$se_sem_tail),
      lambda    = c(res$lambda_sem_wing,    res$lambda_sem_tail),
      se_lambda = c(res$se_lambda_sem_wing, res$se_lambda_sem_tail),
      p.value   = NA_real_,
      Approach  = c("SEM_wing", "SEM_tail")
    )
  }) %>%
    list_rbind(names_to = "Species") %>%
    mutate(LCI95 = estimate - 1.96 * std.error,
           UCI95 = estimate + 1.96 * std.error)

  print(sem_parms)

} # end if(include_sem)


# Plot slope estimates ----------------------------------------------------
approach_labs <- c(
  "Ratio"     = "Wing / Mass",
  "Ratio2"    = "Wing² / Mass",
  "Sli_est"   = "SLI estimated",
  "Sli_iso"   = "SLI isometry",
  "Resid_ols" = "OLS residuals",
  "Ryding"    = "Mass as covariate",
  "SEM_wing"  = "SEM (latent Size)"
)

plot_data <- parms_df
if (include_sem) plot_data <- bind_rows(plot_data, sem_parms %>% filter(Approach == "SEM_wing"))

plot_data %>% filter(term == "B.Lat") %>%
  mutate(Approach = factor(Approach),
         Approach = fct_reorder(.f = Approach, .x = estimate, .fun = mean, .desc = TRUE)) %>%
  ggplot(aes(x = Species, y = estimate, color = Approach,
             group = interaction(Species, Approach))) +
  geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                alpha = .8, width = 0,
                position = position_dodge(width = 0.75)) +
  geom_point(size = 2, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = expression(beta[Lat] ~ "on wing shape")) +
  scale_color_hue(labels = approach_labs)

ggsave("Figures/Nightjar_shape.png", bg = "white")


# SEM diagnostics: factor loadings ----------------------------------------
# Only when include_sem = TRUE.

if (include_sem) {

  # lambda = loading of each appendage on latent Size.
  # low lambda  → appendage barely tracks body size; latent factor poorly
  #               constrained → high SE. The collider problem barely exists.
  # high lambda → appendage dominated by size variance; little independent
  #               variance left for the direct B.Lat path → also high SE.
  # optimal     → intermediate lambda gives best precision.

  sem_diag <- sem_parms %>%
    filter(!is.na(lambda)) %>%
    mutate(
      appendage = if_else(Approach == "SEM_wing", "Wing", "Tail"),
      ci_width  = UCI95 - LCI95
    ) %>%
    left_join(Spp_metadata, by = "Species")

  # 1. Lambda vs SE of B.Lat coefficient (feasibility diagnostic)
  p_lambda_se <- ggplot(sem_diag, aes(x = lambda, y = std.error, color = appendage)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_point(size = 3) +
    geom_text(aes(label = Species), hjust = -0.15, size = 3) +
    scale_color_manual(values = c("Wing" = "steelblue", "Tail" = "coral3")) +
    annotate("text", x = -Inf, y = 0.52, hjust = -0.1, size = 3,
             label = "SE = 0.5", color = "grey40") +
    labs(x = expression(lambda ~ "(factor loading on latent Size)"),
         y = "SE of B.Lat coefficient",
         color = NULL)

  # 2. Lambda vs coefficient estimate with 95% CI
  p_lambda_coef <- ggplot(sem_diag, aes(x = lambda, y = estimate, color = appendage)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95), width = 0, alpha = 0.35) +
    geom_point(size = 3) +
    geom_text(aes(label = Species), hjust = -0.15, size = 3) +
    scale_color_manual(values = c("Wing" = "steelblue", "Tail" = "coral3")) +
    labs(x = expression(lambda ~ "(factor loading on latent Size)"),
         y = "B.Lat coefficient\n(negative = Allen's rule)",
         color = NULL)

  # 3. Lambda vs its own SE (identifies imprecisely estimated loadings)
  p_lambda_uncertainty <- ggplot(sem_diag, aes(x = lambda, y = se_lambda, color = appendage)) +
    geom_point(size = 3) +
    geom_text(aes(label = Species), hjust = -0.15, size = 3) +
    scale_color_manual(values = c("Wing" = "steelblue", "Tail" = "coral3")) +
    labs(x = expression(lambda ~ "(factor loading on latent Size)"),
         y = expression("SE of" ~ lambda),
         color = NULL)

  plot_grid(p_lambda_se, p_lambda_coef, p_lambda_uncertainty, nrow = 1)
  ggsave("Figures/Nightjar_sem_diagnostics.png", bg = "white", width = 14, height = 5)

  # SE(lambda) vs SE(temperature coefficient): just-identification signature.
  # Both SEs reflect the same likelihood flatness — the model cannot cleanly
  # separate the indirect (Size-mediated) from the direct (Allen) B.Lat path.
  sem_diag %>%
    ggplot(aes(x = se_lambda, y = std.error, color = appendage)) +
    geom_point(size = 3) +
    geom_text(aes(label = Species), hjust = -0.1, size = 3) +
    scale_color_manual(values = c("Wing" = "steelblue", "Tail" = "coral3")) +
    labs(x = expression("SE of" ~ lambda),
         y = "SE of B.Lat coefficient",
         color = NULL)

  # Sample size vs SE(lambda)
  sem_diag %>%
    ggplot(aes(x = num_obs, y = se_lambda, color = appendage)) +
    geom_point(size = 3) +
    geom_text(aes(label = Species), hjust = -0.1, size = 3) +
    scale_color_manual(values = c("Wing" = "steelblue", "Tail" = "coral3")) +
    geom_vline(xintercept = 150, linetype = "dashed") +
    labs(x = "Sample size (Wing + Mass + Tail complete)",
         y = expression("SE of" ~ lambda),
         color = NULL)

  # OLS allometric slope vs SE(lambda): weak allometric coupling = poorly
  # identified latent factor
  sem_diag %>%
    ggplot(aes(x = b_ols, y = se_lambda, color = appendage)) +
    geom_point(size = 3) +
    geom_text(aes(label = Species), hjust = -0.1, size = 3) +
    scale_color_manual(values = c("Wing" = "steelblue", "Tail" = "coral3")) +
    labs(x = "OLS allometric slope (log Wing ~ log Mass)",
         y = expression("SE of" ~ lambda),
         color = NULL)

  # Factor loading summary table (sorted by SE descending — worst cases first)
  sem_diag %>%
    dplyr::select(Species, appendage, lambda, se_lambda, estimate, std.error,
                  LCI95, UCI95, num_obs, b_ols, cor_mw) %>%
    arrange(desc(std.error))

} # end if(include_sem)
