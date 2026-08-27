## Analysis -- Shapeshifting in nightjars ##


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(smatr)
library(cowplot)
library(broom)
library(ggpmisc)
select <- dplyr::select
ggplot2::theme_set(theme_cowplot())

library(sliR)   # SLI + simulation functions; see github.com/LezzGitIt/sliR
source("Scripts/00_Key_allometry_fns.R")
nj_raw <- read.csv("Data/Capri_BA_compare03.29.26.csv")

# Control parameters --------------------------------------------------------
control_age_sex <- TRUE   # TRUE: include age / sex as covariates where they significantly affect body size
just_am         <- FALSE  # restrict entire analysis to Adult Males only
# (Age == "Adult" & Sex == "M"); incompatible with control_age_sex
p_age_sex       <- 0.10   # p-value threshold: age / sex covariate inclusion
min_n_age_group <- 100    # min individuals per age group to include age control
min_n_sex_group <- 100    # min individuals per sex group to include sex control
cor_min         <- 0.3    # min Pearson r (mass ~ wing) within group for allometric correlation filter
cor_p_max       <- 0.05   # max p-value for mass ~ wing within group
pos_allom       <- TRUE   # retain only species with cor_mw >= cor_min & p_mw < cor_p_max for SLI-estimated

if (just_am && control_age_sex) {
  stop("`just_am = TRUE` is incompatible with `control_age_sex`. Update control flags")
}
# Age and Sex are held to the same per-group sample-size bar everywhere they're used as a
# single-covariate grouping (Ryding/OLS-residual's covariate decision, and
# calc_sli_hierarchical()'s Level-2 marginal check) -- min_n_age_group is passed explicitly
# as that function's n_min_marginal, so this assertion is what actually keeps them tied
# together rather than the two just happening to be set to the same number.
stopifnot("min_n_age_group and min_n_sex_group must match -- both feed calc_sli_hierarchical()'s single n_min_marginal argument" =
            min_n_age_group == min_n_sex_group)

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

## Pearson r (not a regression slope, so it's on the same [-1,1] scale as cor_min
## below) between Mass and Wing, matching Weeks_2020_ral.R / Atlantic_birds_shape.R.
Cors_tbl <- map(nj_df_analysis_l, \(df){
  ct <- cor.test(df$Mass, df$Wing, use = "complete.obs")
  tibble(cor_mw = as.numeric(ct$estimate), p_mw = ct$p.value)
}) %>% list_rbind(names_to = "Species")

Spp_metadata <- left_join(Spp_metadata, Cors_tbl, by = "Species")
Spp_metadata

# Species-level allometric-correlation filter (same as Weeks / Atlantic): gates
# whether SLI-estimated is computed at all for a species (sli_estimated = NA
# otherwise), since SMA is unreliable when the mass-wing correlation is weak
# (Smith 2009).
Spp_metadata <- Spp_metadata %>%
  mutate(Keep = ifelse(cor_mw >= cor_min & p_mw < cor_p_max, "Include", "Exclude"))

Spp_keep_vec <- Spp_metadata %>%
  filter(Keep == "Include") %>%
  pull(Species)
Spp_keep_vec

if (pos_allom) message(length(Spp_keep_vec), " / ", nrow(Spp_metadata), " species pass allometric filter (sli_est will be NA for the rest)")

# Species-level mass~wing correlation (exported for manuscript) -----------
# Used in Discussion, Limitations of the SLI approaches, to report what
# fraction of species per study fall below the cor_min allometric-reliability
# threshold.
Allometric_cor_nj <- Spp_metadata %>%
  dplyr::select(species = Species, cor_mw, p_mw, Keep) %>%
  mutate(Study = "Nightjar")
write_csv(Allometric_cor_nj, "Derived/Csv/Nightjar_allometric_cor.csv")

# Age / Sex significance testing -------------------------------------------
# Test whether Age / Sex significantly affect Mass or Wing, controlling for
# B.Temp -- only species where Age or Sex show a significant effect (p < p_age_sex)
# are considered for age/sex-stratified covariates or per-group SMA slopes below
# (same logic as Weeks_2020_ral.R / Atlantic_birds_shape.R's sig_age_any/sig_sex_any).
Age_tbl <- bind_rows(
  test_group_effect(nj_df_analysis_l, dv = "Mass", iv = "Age", gradient = "B.Temp", min_n_per_group = min_n_age_group),
  test_group_effect(nj_df_analysis_l, dv = "Wing", iv = "Age", gradient = "B.Temp", min_n_per_group = min_n_age_group)
) %>% mutate(sig = p.value < p_age_sex)

Sex_tbl <- bind_rows(
  test_group_effect(nj_df_analysis_l, dv = "Mass", iv = "Sex", gradient = "B.Temp", min_n_per_group = min_n_sex_group),
  test_group_effect(nj_df_analysis_l, dv = "Wing", iv = "Sex", gradient = "B.Temp", min_n_per_group = min_n_sex_group)
) %>% mutate(sig = p.value < p_age_sex)

Age_tbl %>% filter(sig) %>% dplyr::select(species_, dv, estimate, p.value)
Sex_tbl %>% filter(sig) %>% dplyr::select(species_, dv, estimate, p.value)

sig_age_mass <- Age_tbl %>% filter(dv == "Mass", sig) %>% pull(species_)
sig_age_wing <- Age_tbl %>% filter(dv == "Wing", sig) %>% pull(species_)
sig_sex_mass <- Sex_tbl %>% filter(dv == "Mass", sig) %>% pull(species_)
sig_sex_wing <- Sex_tbl %>% filter(dv == "Wing", sig) %>% pull(species_)

# Union: if age/sex affects either morphometric, include in all downstream models
sig_age_any <- unique(c(sig_age_mass, sig_age_wing))
sig_sex_any <- unique(c(sig_sex_mass, sig_sex_wing))

# Nighthawk was not reliably aged (Age = "Unk" for the large majority of records) --
# exclude it from Age-covariate consideration outright, rather than relying on
# min_n_age_group to filter it out incidentally (its small Adult/Young subgroups
# already fall below that threshold today, but this is a data-quality exclusion,
# not a sample-size one, so it shouldn't depend on that threshold staying put).
sig_age_any <- setdiff(sig_age_any, "Nighthawk")

if (!control_age_sex) {
  sig_age_any <- character(0)
  sig_sex_any <- character(0)
}

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

# Per species: no individuals are excluded from df for age/sex-group reliability reasons --
# the full species dataset feeds all six approaches (see Project_notes.md for the
# rationale). Two separate, independently-gated uses of Age/Sex follow:
#   (1) covs: Ryding/OLS-residual's covariate decision -- Age/Sex enters those two
#       regressions' formula only when >=2 of its levels each have >= min_n_*_group
#       individuals, a strict, species-level, all-or-nothing decision. No mass-wing
#       correlation requirement here: a plain factor covariate in an OLS regression doesn't
#       need a strong within-group correlation to be well-estimated -- that requirement is
#       specific to SMA slopes (SMA slope = OLS slope / r, so a weak r makes the slope
#       numerically unstable; Smith 2009), which is why it still gates (2) below. Unknown-
#       coded individuals are excluded only from this decision (not from df); if the
#       covariate is included, they pick up their own "Unk"/"U" factor level in the fitted
#       regression like any other level.
#   (2) covs_sli (below): SLI-estimated's per-individual hierarchical cascade
#       (calc_sli_hierarchical()), gated only by whether the covariate was significant at
#       all (sig_age_any/sig_sex_any) -- looser than covs, since the cascade itself
#       determines per individual whether a cell/marginal group is reliable enough to use
#       (including the r>=cor_min/p<cor_p_max SMA-reliability check), falling back to the
#       species-wide pooled slope otherwise.
nj_df_l2 <- imap(nj_df_analysis_l, \(df, sp) {
  covs <- character(0)

  if (sp %in% sig_age_any) {
    valid_age <- df %>% filter(!is.na(Age) & Age != "Unk") %>%
      count(Age) %>% filter(n >= min_n_age_group) %>% pull(Age)
    if (length(valid_age) >= 2) covs <- c(covs, "Age")
  }
  if (sp %in% sig_sex_any) {
    valid_sex <- df %>% filter(!is.na(Sex) & Sex != "U") %>%
      count(Sex) %>% filter(n >= min_n_sex_group) %>% pull(Sex)
    if (length(valid_sex) >= 2) covs <- c(covs, "Sex")
  }

  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  ols_mod   <- lm(as.formula(paste("log_wing ~ log_mass", covs_str)), data = df)
  sma_mod   <- sma(log_wing ~ log_mass, data = df, method = "SMA")
  est_b_sma <- coef(sma_mod)["slope"]

  df <- df %>%
    mutate(resid_ols = log_wing - predict(ols_mod, newdata = df),
           resid_sma = residuals(sma_mod))

  ## Estimated SLI: per-group SMA slopes when this species has valid covariates, otherwise the species-wide SMA slope. Kept as separate calls because sliR::calc_sli() ignores b_sli whenever control is supplied, so passing both would silently discard one of them.
  df_iso <- df %>% sliR::calc_sli(b_sli = 0.33, Append = Wing, rename_col = "sli_isometry")

  ## SLI-estimated: calc_sli_hierarchical() now handles the species-level reliability gate
  ## internally (returns NA for every individual if the species-wide correlation itself is
  ## unreliable) and falls back to the plain pooled slope when no covariate was significant
  ## (covs_sli empty), so one unconditional call covers every species -- matches Spp_keep_vec
  ## exactly, since both check the same species-wide mass~wing correlation.
  covs_sli <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
  calc_sli_hierarchical(df_iso, Append = Wing, covariates = covs_sli, rename_col = "sli_estimated",
                        n_min_marginal = min_n_age_group)
})

# Per-group allometric correlation (Mass ~ Wing within each Age × Sex combination) --
# Inspect r_mw and p_mw per group; groups with pass = FALSE lack a meaningful
# allometric relationship and should not drive per-group SMA slope estimates.
if (control_age_sex) {
  group_cor_wing <- imap(nj_df_analysis_l, \(df, sp) {
    covs <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
    if (!length(covs)) return(NULL)
    build_group_cor_tbl(df, Append = Wing, Mass = Mass, control = covs) %>%
      mutate(Species = sp, .before = 1)
  }) %>% list_rbind() %>%
    mutate(pass = r_mw >= cor_min & p_mw < cor_p_max)
  print(group_cor_wing)
}

# Per-species per-group SMA slope summary (only for species where Age/Sex is
# actually used as a covariate; see sig_age_any/sig_sex_any above).
if (control_age_sex) {
  sli_slopes_tbl <- imap(nj_df_analysis_l, \(df, sp) {
    covs <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
    if (!length(covs)) return(NULL)
    sliR::build_sli_slopes_tbl(df, Append = Wing, control = covs)
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
    tibble(Metric = "Ratio2", n = nrow(df), r = as.numeric(ct2$estimate), p_value = ct2$p.value),
    build_sli_mass_cor_tbl(df, Mass = Mass)
  )
}) %>% list_rbind(names_to = "Species") %>%
  mutate(Study = "Nightjar", species = Species)

write_csv(Ratio_mass_cor_nj, "Derived/Csv/Nightjar_mass_cor.csv")

# Run models & extract parms ----------------------------------------------
# Approach-specific covariate rules (matching Weeks_2020_ral.R / Atlantic_birds_shape.R):
#   Ratio/Ratio2/Sli_iso: no Age/Sex by design
#   Resid_ols: Age/Sex cleaned in first model (nj_df_l2); no additional covariates
#   Ryding: include Age/Sex in combined model (B.Temp conditional on both), only for
#     species where Age/Sex showed a significant effect (sig_age_any/sig_sex_any)
#   Sli_est: per-group SMA slopes handled upstream via sliR::calc_sli(control = covs);
#     NA for species failing the allometric-correlation filter (Spp_keep_vec)
parms_df <- imap(nj_df_l3, \(df, sp) {
  covs     <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
  covs     <- covs[vapply(covs, \(v) length(unique(na.omit(df[[v]]))) >= 2, logical(1))]
  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  mod_resid_ols   <- lm(resid_ols    ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Resid_ols")
  mod_coef_ols    <- lm(as.formula(paste("Wing ~ Mass + B.Temp", covs_str)), data = df) %>%
    tidy() %>% mutate(Approach = "Ryding")
  mod_coef_ratio  <- lm(wing_mass    ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass   ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Ratio2")
  mod_sli_iso     <- lm(sli_isometry ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Sli_iso")
  mod_sli_est <- if (any(!is.na(df$sli_estimated))) {
    lm(sli_estimated ~ B.Temp, data = df) %>% tidy() %>% mutate(Approach = "Sli_est")
  } else tibble()

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

# Strict rank-order check (no family averaging): A/S > A2/S > SLI-isometry > multiple regression
# for Bergmann's species (reversed for Inverse Bergmann's), following each method's implicit
# scaling exponent directly (A/S: beta=1, A2/S: beta=0.5, SLI-isometry: beta=0.33) rather than
# family-averaged coefficients. SLI-estimated is excluded because its exponent is each species'
# own empirical SMA slope, not a fixed constant, so it has no fixed position in this chain;
# OLS-residuals is excluded because it isn't part of the proposed ordering either.
rank_strict_nj <- parms_df %>%
  filter(term == "B.Temp", Approach %in% c("Ratio", "Ratio2", "Sli_iso", "Ryding")) %>%
  dplyr::select(Species, Approach, estimate) %>%
  pivot_wider(names_from = Approach, values_from = estimate) %>%
  left_join(Direction_nj, by = "Species") %>%
  mutate(
    rank_strict_consistent = case_when(
      Direction == "Bergmann's"         ~ Ratio > Ratio2 & Ratio2 > Sli_iso & Sli_iso > Ryding,
      Direction == "Inverse Bergmann's" ~ Ratio < Ratio2 & Ratio2 < Sli_iso & Sli_iso < Ryding,
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(Species, rank_strict_consistent)

# CSV export ---------------------------------------------------------------
nj_parms_out <- parms_df %>%
  filter(term == "B.Temp") %>%
  left_join(Direction_nj, by = "Species") %>%
  left_join(rank_nj, by = "Species") %>%
  left_join(rank_strict_nj, by = "Species") %>%
  mutate(
    Study    = "Nightjar",
    species_ = str_replace_all(Species, " ", "_"),
    species  = Species
  )
write_csv(nj_parms_out, "Derived/Csv/Nightjar_parms.csv")

