## Analysis -- Shapeshifting in nightjars ##


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(smatr)
library(cowplot)
library(broom)
library(ggpmisc)
select <- dplyr::select
ggplot2::theme_set(theme_cowplot())

source("Scripts/Key_allometry_fns.R")
nj_raw <- read.csv("Data/Capri_BA_compare03.29.26.csv")

# Control parameters --------------------------------------------------------
control_age_sex_6mod <- TRUE  # include Age/Sex control in the six OLS-based approaches
# Nighthawk: Sex only; Nightjar/Whip-poor-will: Age + Sex
just_am              <- FALSE   # restrict entire analysis to Adult Males only
# (Age == "Adult" & Sex == "M"); incompatible with control_age_sex_6mod

if (just_am && control_age_sex_6mod) {
  stop("`just_am = TRUE` is incompatible with `control_age_sex_6mod`. Update control flags")
}

# Download temperature data -------------------------------------------------
# WorldClim Annual Max Temperature, °C
# Cached to Data/Nightjar_temp.rds after first run; subsequent runs skip the
# WorldClim download entirely. Delete the cache file to force a refresh.
nj_temp_cache <- "Data/Nightjar_temp.rds"
if (file.exists(nj_temp_cache)) {
  nj_raw <- readRDS(nj_temp_cache)
} else {
  library(geodata)
  library(terra)
  tmax <- worldclim_global(var = "tmax", res = 2.5, path = "Data/")
  mean_tmax <- mean(tmax[[5:9]])
  coords <- cbind(nj_raw$B.Long, nj_raw$B.Lat)
  plot(mean_tmax)
  nj_raw$B.Temp <- terra::extract(mean_tmax, coords)[, 1]  # single-layer raster → 1-col output
  saveRDS(nj_raw, nj_temp_cache)
}

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
# NOTE on sign convention: B.Temp is annual mean temperature (°C, WorldClim
# BIO1) at the banding location. A POSITIVE B.Temp coefficient indicates
# Allen's rule (longer appendages at warmer sites). This matches the simulation
# convention where higher Temp = warmer.

nj_df <- nj_raw %>%
  dplyr::select(Species,
                Wing = Wing.comb, Mass = Mass.comb, Tail = Tail.comb,
                B.Temp, Age, Sex) %>%
  drop_na(Wing, Mass) %>%
  mutate(log_wing   = log(Wing),
         log_mass   = log(Mass),
         log_tail   = log(Tail),
         wing_mass  = log_wing - log_mass,
         wing2_mass = 2*log_wing - log_mass)

if (just_am) nj_df <- nj_df %>% filter(Age == "Adult" & Sex == "M")

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
# Wing + Mass complete cases — analysis dataset for all six OLS/SMA/SLI/ratio methods
nj_df_l <- nj_df %>% group_split(Species)
names(nj_df_l) <- c("Nighthawk", "Nightjar", "Whip-poor-will")

nj_df_analysis_l <- nj_df_l

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
  lm(log_wing ~ log_mass + B.Temp, data = df)
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

# Per-species covariate helper for the six OLS-based approaches.
# Nighthawk: Sex only (Age is mostly "Unk" for this species).
# Nightjar / Whip-poor-will: Age + Sex.
# Returns character(0) when control_age_sex_6mod = FALSE (no control applied).
get_6mod_covs <- function(sp) {
  if (!control_age_sex_6mod) return(character(0))
  if (sp == "Nighthawk") return("Sex")
  c("Age", "Sex")
}

nj_df_l2 <- imap(nj_df_analysis_l, \(df, sp) {
  covs     <- get_6mod_covs(sp)
  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  # Drop unknown-coded rows early so predict() never encounters new factor levels.
  # For Nighthawk (covs = "Sex"): removes Sex == "U" / NA only.
  # For Nightjar / Whip-poor-will (covs = c("Age","Sex")): also removes Age == "Unk" / NA.
  if (length(covs)) {
    df <- df %>% filter(if_all(all_of(covs), \(x) !is.na(x) & x != "Unk" & x != "U"))
  }

  ols_mod   <- lm(as.formula(paste("log_wing ~ log_mass", covs_str)), data = df)
  sma_mod   <- sma(log_wing ~ log_mass, data = df, method = "SMA")
  est_b_sma <- coef(sma_mod)["slope"]

  df <- df %>%
    mutate(resid_ols = log_wing - predict(ols_mod, newdata = df),
           resid_sma = residuals(sma_mod))

  ## Estimated SLI: per-group SMA slopes when this species has valid covariates, otherwise the species-wide SMA slope. Kept as separate calls because calc_sli() ignores b_sli whenever control is supplied, so passing both would silently discard one of them.
  df_iso <- df %>% calc_sli(b_sli = 0.33, Append = Wing, rename_col = "sli_isometry")

  if (length(covs)) {
    calc_sli(df_iso, Append = Wing, control = covs, rename_col = "sli_estimated")
  } else {
    calc_sli(df_iso, Append = Wing, b_sli = est_b_sma, rename_col = "sli_estimated")
  }
})

# Per-species per-group SMA slope summary (10 rows total when control_age_sex_6mod = TRUE):
# 2 rows for Nighthawk (Sex), 4 for Nightjar (Age x Sex), 4 for Whip-poor-will (Age x Sex).
if (control_age_sex_6mod) {
  sli_slopes_tbl <- imap(nj_df_analysis_l, \(df, sp) {
    covs <- get_6mod_covs(sp)
    if (!length(covs)) return(NULL)
    build_sli_slopes_tbl(df, Append = Wing, control = covs)
  }) %>% list_rbind(names_to = "Species")
  print(sli_slopes_tbl)
}

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

# Ratio-mass correlation (exported for manuscript) --------------------------
# Correlation between the wing/mass ratio (log(A/S) under log_ratio = TRUE)
# and its squared counterpart (log(A^2/S)) with body mass, illustrating the
# confounding of ratio metrics with body size and how squaring the appendage
# term changes it.
Ratio_mass_cor_nj <- imap(nj_df_l3, \(df, sp) {
  ct1 <- cor.test(df$wing_mass,  df$Mass)
  ct2 <- cor.test(df$wing2_mass, df$Mass)
  bind_rows(
    tibble(Metric = "Ratio",  n = nrow(df), r = as.numeric(ct1$estimate), p_value = ct1$p.value),
    tibble(Metric = "Ratio2", n = nrow(df), r = as.numeric(ct2$estimate), p_value = ct2$p.value)
  )
}) %>% list_rbind(names_to = "Species") %>%
  mutate(Study = "Nightjar", species = Species)

write_csv(Ratio_mass_cor_nj, "Derived/Csv/Nightjar_ratio_mass_cor.csv")

# Run models & extract parms ----------------------------------------------
parms_df <- imap(nj_df_l3, \(df, sp) {
  # df is already filtered (Unk/U rows dropped in nj_df_l2) so no further subsetting needed.
  covs     <- get_6mod_covs(sp)
  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  # Resid_ols: residuals already age/sex-cleaned via first OLS model in nj_df_l2.
  mod_resid_ols   <- lm(resid_ols    ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Resid_ols")
  # Ryding: age/sex as covariates in the combined model (B.Temp conditional on both).
  mod_coef_ols    <- lm(as.formula(paste("Wing ~ Mass + B.Temp", covs_str)), data = df) %>%
    tidy() %>% mutate(Approach = "Ryding")
  # Ratio/Ratio2/Sli_iso: no age/sex control by design.
  mod_coef_ratio  <- lm(wing_mass    ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass   ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Ratio2")
  mod_sli_iso     <- lm(sli_isometry ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Sli_iso")
  # Sli_est: per-group SMA slopes (from calc_sli) already handle age/sex variation;
  # no covariates in the final regression. NA rows dropped automatically.
  mod_sli_est     <- lm(sli_estimated ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Sli_est")

  bind_rows(mod_coef_ratio, mod_coef_ratio2, mod_coef_ols, mod_resid_ols, mod_sli_est, mod_sli_iso)
}) %>% list_rbind(names_to = "Species") %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

# Plot slope estimates ----------------------------------------------------
approach_labs <- c(
  "Ratio"     = "Wing / Mass",
  "Ratio2"    = "Wing² / Mass",
  "Sli_est"   = "SLI estimated",
  "Sli_iso"   = "SLI isometry",
  "Resid_ols" = "OLS residuals",
  "Ryding"    = "Mass as covariate"
)

parms_df %>% filter(term == "B.Temp") %>%
  mutate(Approach = factor(Approach, levels = c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso"))) %>%
  ggplot(aes(x = Species, y = estimate, color = Approach,
             group = interaction(Species, Approach))) +
  geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                alpha = .8, width = 0,
                position = position_dodge(width = 0.75)) +
  geom_point(size = 2, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = expression(beta[T] ~ "on wing shape")) +
  scale_color_hue(labels = approach_labs)

ggsave("Figures/Nightjar_shape.png", bg = "white")


# Direction classification -------------------------------------------------
dir_mods_nj <- imap(nj_df_analysis_l, \(df, sp) {
  bind_rows(
    lm(Mass ~ B.Temp, data = df) %>% tidy() %>% filter(term == "B.Temp") %>%
      mutate(Species = sp, dv = "mass"),
    lm(Wing ~ B.Temp, data = df) %>% tidy() %>% filter(term == "B.Temp") %>%
      mutate(Species = sp, dv = "wing")
  )
}) %>% list_rbind()

Direction_nj <- classify_direction(dir_mods_nj, species_col = "Species",
                                   mass_dv = "mass", wing_dv = "wing")

# Rank consistency (same logic as Weeks / Atlantic) -----------------------
rank_nj <- parms_df %>%
  filter(term == "B.Temp",
         Approach %in% c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso")) %>%
  mutate(approach_group = case_when(
    Approach %in% c("Ratio", "Ratio2")     ~ "avg_ratio",
    Approach %in% c("Sli_est", "Sli_iso")  ~ "avg_sli",
    Approach %in% c("Ryding", "Resid_ols") ~ "avg_ols"
  )) %>%
  group_by(Species, approach_group) %>%
  summarise(avg_coef = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = approach_group, values_from = avg_coef) %>%
  left_join(Direction_nj, by = "Species") %>%
  mutate(
    rank_consistent = case_when(
      Direction == "Bergmann's"         ~ avg_ratio > avg_sli & avg_sli > avg_ols,
      Direction == "Inverse Bergmann's" ~ avg_ratio < avg_sli & avg_sli < avg_ols,
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(Species, rank_consistent)

# CSV export ---------------------------------------------------------------
nj_parms_out <- parms_df %>%
  filter(term == "B.Temp") %>%
  left_join(Direction_nj, by = "Species") %>%
  left_join(rank_nj, by = "Species") %>%
  mutate(
    Study    = "Nightjar",
    species_ = str_replace_all(Species, " ", "_"),
    species  = Species
  )
write_csv(nj_parms_out, "Derived/Csv/Nightjar_parms.csv")

