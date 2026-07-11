## Examining changes in wing shape in the Brazilian Atlantic Forest ##

# Libraries & data -------------------------------------------------
library(tidyverse)
library(broom)
library(janitor)
library(smatr)
library(cowplot)
library(patchwork)
library(naniar)
library(sliR)   # SLI + simulation functions; see github.com/LezzGitIt/sliR
source("Scripts/Key_allometry_fns.R")
ggplot2::theme_set(theme_cowplot())

Atlantic_birds <- read_csv("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Ecology/Atlantic_bird_traits/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv")

# Run parameters ----------------------------------------------------------
min_n_obs       <- 150    # minimum observations per species to be included
year_cutoff     <- 1990   # exclude records from this year and before
p_bergmann      <- 0.05   # p-value threshold: Bergmann / Allen response classification
p_age_sex       <- 0.10   # p-value threshold: age / sex covariate inclusion
pos_allom       <- TRUE  # retain only species with cor_mw > cor_min & p_mw < cor_p_max
filter_bergmann <- FALSE  # TRUE: retain only species with a significant Bergmann / Allen response
control_age_sex <- TRUE   # TRUE: include age / sex as covariates where they significantly affect body size
min_n_age_group <- 100    # min individuals per age group (juvenile/adult) to include age control
min_n_sex_group <- 100    # min individuals per sex group (female/male) to include sex control
cor_min         <- 0.3    # min Pearson r (mass ~ wing) within group for allometric correlation filter
cor_p_max       <- 0.05   # max p-value for mass ~ wing within group

# Format ------------------------------------------------------------------
# Combine several metrics of wing, filter out observations with no wing or mass
Atlantic_birds2 <- Atlantic_birds %>% 
  clean_names() %>% 
  rename(mass = body_mass_g, body_length = body_length_mm, B.Tavg = annual_mean_temperature, species_epithet = species, latitude = latitude_decimal_degrees) %>%
  mutate(species_ = str_replace(binomial, " ", "_"), 
         wing = coalesce(
           wing_length_mm, wing_length_left_mm, wing_length_right_mm
           ),
         tarsus = coalesce(tarsus_length_mm, tarsus_length_right_mm, tarsus_length_left_mm),
         tail = coalesce(tail_length_mm, tail_length_right_mm, tail_length_left_mm)
         )

# Format predictor & responsible variables 
Atlantic_birds3 <- Atlantic_birds2 %>% 
  mutate(log_wing = log(wing), 
         log_mass = log(mass),
         log_tarsus = log(tarsus),
         log_tail = log(tail),
         wing_mass = log_wing - log_mass,
         wing2_mass = 2*log_wing - log_mass)

# Morphology -------------------------------------------------------------
# Inspect possible measurements: Body length, wing, mass, tarsus, tail
Atlantic_birds3 %>% 
  dplyr::select(
   tarsus, wing, tail, mass, body_length, wingspan_mm, head_length_total_mm, bill_width_mm, bill_depth_mm, age, sex # starts_with("Bill_")
  ) %>%
  naniar::gg_miss_var(show_pct = TRUE)

# Filter NA morphologies
Atlantic_birds4 <- Atlantic_birds3 %>%
  filter(!is.na(mass) & !is.na(wing) & !is.na(tarsus)) %>% 
  filter(year > year_cutoff & !is.na(B.Tavg)) %>% #sex == "Male" & age == "Adult" &
  dplyr::select(species_, wing, tarsus, tail, mass, log_wing, log_tarsus, log_tail, log_mass, wing_mass, wing2_mass, B.Tavg, latitude, locality, year, age, sex)
table(Atlantic_birds4$year)

# 780 species total
Atlantic_birds3 %>% pull(species_) %>% 
  unique()

# Filter species based on latitudinal range and number of observations
Spp_summary <- Atlantic_birds4 %>%
  summarize(min_lat = min(latitude, na.rm = T),
            max_lat = max(latitude, na.rm = T),
            range_lat = min_lat - max_lat,
            num_obs = n(),
            num_locality = length(unique(locality)),
            .by = species_)
Spp_include <- Spp_summary %>%
  filter(num_obs > min_n_obs) %>%
  arrange(desc(num_obs))

# Create list 
Atl_birds_l <- Atlantic_birds4 %>% 
  right_join(Spp_include) %>% 
  group_split(species_)

# Name list
names(Atl_birds_l) <- map_chr(Atl_birds_l, ~unique(pull(.x, species_)))
length(Atl_birds_l)

# Age & sex --------------------------------------------------------------
Atlantic_birds4 %>% tabyl(age)
Atlantic_birds4 %>% tabyl(sex)

# Test effect of a grouping variable (iv = "age" or "sex") on a morphometric DV.
# Filters Unknown/NA for iv, then drops individual groups below min_n_per_group.
# If fewer than 2 valid groups remain, returns NULL for that species.
test_group_effect <- function(df_list, dv, iv, min_n_per_group = 0) {
  map(df_list, \(df) {
    df_filt <- df %>% filter(.data[[iv]] != "Unknown" & !is.na(.data[[iv]]))
    if (nrow(df_filt) < 10) return(NULL)

    valid_grps <- df_filt %>%
      count(.data[[iv]]) %>%
      filter(n >= min_n_per_group) %>%
      pull(.data[[iv]])
    if (length(valid_grps) < 2) return(NULL)
    df_filt <- df_filt %>% filter(.data[[iv]] %in% valid_grps)

    n_counts <- df_filt %>%
      count(.data[[iv]], name = "n") %>%
      mutate(lbl = paste0("n_", tolower(.data[[iv]]))) %>%
      dplyr::select(lbl, n) %>%
      pivot_wider(names_from = lbl, values_from = n)

    fmla <- as.formula(paste(dv, "~ B.Tavg +", iv))
    tidy(lm(fmla, data = df_filt)) %>%
      filter(str_starts(term, iv)) %>%
      mutate(dv = dv) %>%
      bind_cols(n_counts)
  }) %>%
    list_rbind(names_to = "species_")
}

# Age and sex effects on mass, wing, and tarsus
# Groups below min_n_age_group / min_n_sex_group are dropped inside the function;
# species with fewer than 2 valid groups are excluded entirely.
Age_tbl <- bind_rows(
  test_group_effect(Atl_birds_l, dv = "mass",   iv = "age", min_n_per_group = min_n_age_group),
  test_group_effect(Atl_birds_l, dv = "wing",   iv = "age", min_n_per_group = min_n_age_group),
  test_group_effect(Atl_birds_l, dv = "tarsus", iv = "age", min_n_per_group = min_n_age_group)
) %>% mutate(sig = p.value < p_age_sex)

Sex_tbl <- bind_rows(
  test_group_effect(Atl_birds_l, dv = "mass",   iv = "sex", min_n_per_group = min_n_sex_group),
  test_group_effect(Atl_birds_l, dv = "wing",   iv = "sex", min_n_per_group = min_n_sex_group),
  test_group_effect(Atl_birds_l, dv = "tarsus", iv = "sex", min_n_per_group = min_n_sex_group)
) %>% mutate(sig = p.value < p_age_sex)

# Inspect: which species × DV combinations are significant?
Age_tbl %>% filter(sig) %>% dplyr::select(species_, dv, estimate, p.value, n_juvenile, n_adult)
Sex_tbl %>% filter(sig) %>% dplyr::select(species_, dv, estimate, p.value, n_female, n_male)

# Species with significant age / sex effects by DV (used in Berg's models and downstream)
sig_age_mass   <- Age_tbl %>% filter(dv == "mass",   sig) %>% pull(species_)
sig_sex_mass   <- Sex_tbl %>% filter(dv == "mass",   sig) %>% pull(species_)
sig_age_wing   <- Age_tbl %>% filter(dv == "wing",   sig) %>% pull(species_)
sig_sex_wing   <- Sex_tbl %>% filter(dv == "wing",   sig) %>% pull(species_)
sig_age_tarsus <- Age_tbl %>% filter(dv == "tarsus", sig) %>% pull(species_)
sig_sex_tarsus <- Sex_tbl %>% filter(dv == "tarsus", sig) %>% pull(species_)

# Union across all DVs: if age/sex affects any indicator, include it in shape models
sig_age_any <- unique(c(sig_age_mass, sig_age_wing, sig_age_tarsus))
sig_sex_any <- unique(c(sig_sex_mass, sig_sex_wing, sig_sex_tarsus))

# If control_age_sex = FALSE, zero out — no covariate adjustment flows through downstream
if (!control_age_sex) {
  sig_age_any <- character(0)
  sig_sex_any <- character(0)
}

# Berg's relationship -------------------------------------------
## For each species, age/sex included as covariates where they significantly affect the DV; those species' data are also filtered to non-Unknown observations for that variable.

berg_model <- function(df, dv, sig_age_spp, sig_sex_spp) {
  sp   <- unique(df$species_)
  df_  <- df
  covs <- character(0)
  if (sp %in% sig_age_spp) {
    df_  <- df_ %>% filter(age != "Unknown" & !is.na(age))
    covs <- c(covs, "age")
  }
  if (sp %in% sig_sex_spp) {
    df_  <- df_ %>% filter(sex != "Unknown" & !is.na(sex))
    covs <- c(covs, "sex")
  }
  rhs  <- paste(c("B.Tavg", covs), collapse = " + ")
  tidy(lm(as.formula(paste(dv, "~", rhs)), data = df_)) %>%
    filter(term == "B.Tavg") %>%
    mutate(dv = dv)
}

mass_mod <- map(Atl_birds_l, berg_model, dv = "mass",
                sig_age_spp = sig_age_mass, sig_sex_spp = sig_sex_mass) %>%
  list_rbind(names_to = "species_")

wing_mod <- map(Atl_birds_l, berg_model, dv = "wing",
                sig_age_spp = sig_age_wing, sig_sex_spp = sig_sex_wing) %>%
  list_rbind(names_to = "species_")

Bergs <- bind_rows(mass_mod, wing_mod)

# Bergs - relationship with temp ----------------------------------
# Select species with a (marginally) significant relationship between temperature and either wing or mass 
# NOTE:: Generally effects are small compared to nightjar dataset 
Spp_keep <- Bergs %>% 
  mutate(Bergs = ifelse(estimate < 0, "Y", "N")) %>%
  #filter(p.value < .5) %>% 
  distinct(species_, Bergs, dv, p.value) %>% 
  arrange(species_)

# Classify the direction of the effect and which trait is significant
Spp_keep2 <- Spp_keep %>%
  mutate(sig = p.value < p_bergmann) %>%
  group_by(species_) %>%
  summarise(
    n_sig     = sum(sig),
    mass_dir  = Bergs[dv == "mass"],   # "Y" = neg estimate (smaller at warmer) = Bergmann's for mass
    wing_dir  = Bergs[dv == "wing"],   # "Y" = neg estimate (shorter at warmer), "N" = longer = Allen's
    mass_sig  = sig[dv == "mass"],
    wing_sig  = sig[dv == "wing"],
    Sig_trait = ifelse(n_sig == 0, "Neither", ifelse(n_sig == 2, "both", dv[sig])),
    .groups   = "drop"
  ) %>%
  mutate(
    sole_dir  = case_when(mass_sig & !wing_sig ~ mass_dir,
                          !mass_sig & wing_sig ~ wing_dir),
    Direction = case_when(
      n_sig == 0                                          ~ "Stable",
      n_sig == 1 & sole_dir == "Y"                       ~ "Bergmann's",
      n_sig == 1 & sole_dir == "N"                       ~ "Inverse Bergmann's",
      n_sig == 2 & mass_dir == "Y" & wing_dir == "Y"     ~ "Bergmann's",
      n_sig == 2 & mass_dir == "N" & wing_dir == "N"     ~ "Inverse Bergmann's",
      n_sig == 2 & mass_dir == "Y" & wing_dir == "N"     ~ "Mixed - Wingier",
      n_sig == 2 & mass_dir == "N" & wing_dir == "Y"     ~ "Mixed - Fatter",
      TRUE ~ "Check"
    )
  )
Spp_keep2 %>% tabyl(Direction)

# If filter_bergmann = TRUE, drop species with no significant Bergmann / Allen response
if (filter_bergmann) Spp_keep2 <- Spp_keep2 %>% filter(Direction != "Stable")

Atlantic_birds5 <- Atlantic_birds4 %>%
  filter(species_ %in% unique(Spp_keep2$species_))

# Calculate correlations, and SMA and OLS slopes
Slopes_tbl <- Atlantic_birds5 %>%
  summarize(b_sma = sd(log_wing, na.rm = T) / sd(log_mass, na.rm = T),
            b_ols = cor(wing, mass) * b_sma,
            .by = species_)

# Examine species remaining
Spp_metadata <- Spp_keep2 %>% 
  left_join(Spp_include) %>% 
  left_join(Slopes_tbl)
Spp_metadata

# Visualize regression approaches -----------------------------------------
## Temp - mass relationship, very small effects
Atlantic_birds5 %>%
  ggplot(aes(x = B.Tavg, y = log_mass, color = species_)) +
  geom_point() +
  geom_smooth(method = "lm") +
  guides(color = "none")

# Compare SMA vs OLS line-fitting methods for allometry
## NOTE: In general that OLS slopes are near 0
Atlantic_birds5 %>%
  ggplot(aes(x = log_mass, y = log_wing, color = species_)) + 
  geom_point(alpha = .2) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, alpha = .3) +
  ggpmisc::stat_ma_line(method = "SMA", se = FALSE, alpha = .3) + 
  #facet_wrap(~species_) + 
  guides(color = "none") 

# Allometric correlation --------------------------------------------------
# Some authors note that SMA should only be used if there is a correlation between the two morphological variables. This makes sense because SMA does not distinguish between scatter (correlation) and the functional relationship between X and Y. See Smith (2009) paper. 

# Create list 
Atl_birds_l2 <- Atlantic_birds5 %>% group_split(species_)

# Name list
names(Atl_birds_l2) <- map_chr(Atl_birds_l2, ~unique(pull(.x, species_)))
length(Atl_birds_l2)

Cors_tbl <- map(Atl_birds_l2, \(df) {
  ct <- cor.test(df$mass, df$wing, use = "complete.obs")
  tibble(cor_mw = as.numeric(ct$estimate), p_mw = ct$p.value)
}) %>% list_rbind(names_to = "species_")
Cors_tbl %>% filter(cor_mw < 0)

Spp_metadata2 <- Spp_metadata %>% left_join(Cors_tbl) %>% 
  mutate(Keep = ifelse(cor_mw > cor_min & p_mw < cor_p_max, "Include", "Exclude"))

Spp_keep_vec <- Spp_metadata2 %>%
  filter(Keep == "Include") %>%
  pull(species_)
Spp_keep_vec

# Atl_birds_l2 is NOT filtered to Spp_keep_vec here; sli_est is set NA for non-passing species in Atl_birds_l3
if (pos_allom) message(length(Spp_keep_vec), " / ", length(Atl_birds_l2), " species pass allometric filter (sli_est will be NA for the rest)")

# Per-group allometric correlation (mass ~ wing within each age × sex combination) -----
# Inspect r_mw and p_mw per group; groups with pass = FALSE lack a meaningful
# allometric relationship and should not drive per-group SMA slope estimates.
if (control_age_sex) {
  group_cor_wing <- imap(Atl_birds_l2, \(df, sp) {
    covs <- c(if (sp %in% sig_age_any) "age", if (sp %in% sig_sex_any) "sex")
    if (!length(covs)) return(NULL)
    build_group_cor_tbl(df, Append = wing, Mass = mass, control = covs) %>%
      mutate(species_ = sp, .before = 1)
  }) %>% list_rbind() %>%
    mutate(pass = r_mw > cor_min & p_mw <= cor_p_max)
  print(group_cor_wing)
}

# Create shape metrics ----------------------------------------------------
# Per species: filter to groups meeting min_n thresholds; include age/sex in
# covariate list only if ≥ 2 valid groups remain after filtering.
# Residuals use predict() so the model and the df are always aligned.
# Sli_est uses per-group SMA slopes via sliR::calc_sli(control = covs); sli_iso uses isometry.
Atl_birds_l3 <- imap(Atl_birds_l2, \(df, sp) {
  covs <- character(0)

  if (sp %in% sig_age_any) {
    valid_age <- df %>%
      filter(!is.na(age) & age != "Unknown") %>%
      count(age) %>% filter(n >= min_n_age_group) %>% pull(age)
    if (length(valid_age) >= 1) df <- df %>% filter(age %in% valid_age)
    if (length(valid_age) >= 2) {
      passing_age <- build_group_cor_tbl(df, Append = wing, Mass = mass, control = "age") %>%
        ungroup() %>% filter(r_mw > cor_min & p_mw <= cor_p_max) %>% pull(age)
      if (length(passing_age) >= 1) df <- df %>% filter(age %in% passing_age)
      if (length(passing_age) >= 2) covs <- c(covs, "age")
    }
  }
  if (sp %in% sig_sex_any) {
    valid_sex <- df %>%
      filter(!is.na(sex) & sex != "Unknown") %>%
      count(sex) %>% filter(n >= min_n_sex_group) %>% pull(sex)
    if (length(valid_sex) >= 1) df <- df %>% filter(sex %in% valid_sex)
    if (length(valid_sex) >= 2) {
      passing_sex <- build_group_cor_tbl(df, Append = wing, Mass = mass, control = "sex") %>%
        ungroup() %>% filter(r_mw > cor_min & p_mw <= cor_p_max) %>% pull(sex)
      if (length(passing_sex) >= 1) df <- df %>% filter(sex %in% passing_sex)
      if (length(passing_sex) >= 2) covs <- c(covs, "sex")
    }
  }

  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  ols_wing     <- lm(as.formula(paste("log_wing   ~ log_mass", covs_str)), data = df)
  ols_tarsus   <- lm(as.formula(paste("log_tarsus ~ log_mass", covs_str)), data = df)
  sma_wing     <- sma(log_wing   ~ log_mass, data = df, method = "SMA")
  sma_tarsus   <- sma(log_tarsus ~ log_mass, data = df, method = "SMA")

  df_res <- df %>%
    mutate(resid_ols        = log_wing   - predict(ols_wing,   newdata = df),
           resid_sma        = residuals(sma_wing),
           resid_ols_tarsus = log_tarsus - predict(ols_tarsus, newdata = df),
           resid_sma_tarsus = residuals(sma_tarsus)) %>%
    sliR::calc_sli(b_sli = 0.33, Append = wing,   Mass = mass, rename_col = "sli_isometry") %>%
    sliR::calc_sli(b_sli = 0.33, Append = tarsus, Mass = mass, rename_col = "sli_tarsus_iso")

  ## Estimated SLI: per-group SMA slopes when this species has valid covariates, otherwise the species-wide SMA slope. Kept as separate branches because sliR::calc_sli() ignores b_sli whenever control is supplied, so passing both would silently discard one of them.
  if (sp %in% Spp_keep_vec) {
    if (length(covs)) {
      df_res %>%
        sliR::calc_sli(Append = wing,   Mass = mass, control = covs, rename_col = "sli_estimated") %>%
        sliR::calc_sli(Append = tarsus, Mass = mass, control = covs, rename_col = "sli_tarsus_est")
    } else {
      df_res %>%
        sliR::calc_sli(Append = wing,   Mass = mass, b_sli = coef(sma_wing)["slope"],   rename_col = "sli_estimated") %>%
        sliR::calc_sli(Append = tarsus, Mass = mass, b_sli = coef(sma_tarsus)["slope"], rename_col = "sli_tarsus_est")
    }
  } else {
    df_res %>% mutate(sli_estimated = NA_real_, sli_tarsus_est = NA_real_)
  }
})

# Per-group SMA slope tables (for inspection — wing and tarsus)
if (control_age_sex) {
  wing_slopes_tbl <- imap(Atl_birds_l2, \(df, sp) {
    covs <- c(if (sp %in% sig_age_any) "age", if (sp %in% sig_sex_any) "sex")
    if (!length(covs)) return(NULL)
    sliR::build_sli_slopes_tbl(df, Append = wing, Mass = mass, control = covs) %>%
      mutate(species_ = sp, .before = 1)
  }) %>% list_rbind()
  print(wing_slopes_tbl)
}

# Scale by species 
Atl_birds_l4 <- map(Atl_birds_l3, \(df){
  df %>% mutate(across(where(is.numeric), scale))
})

# Inspect correlations of shape metric and body size (mass in this case)
map(Atl_birds_l4, \(df){
  df %>% summarize(wm_m = cor(wing_mass, mass),
                   w2m_m = cor(wing2_mass, mass),
                   resid_m = cor(resid_sma, mass))
})

# Ratio-mass correlation (exported for manuscript) --------------------------
# Correlation between the wing/mass ratio (log(A/S) under log_ratio = TRUE)
# and its squared counterpart (log(A^2/S)) with body mass, illustrating the
# confounding of ratio metrics with body size and how squaring the appendage
# term changes it.
Ratio_mass_cor <- imap(Atl_birds_l4, \(df, sp) {
  ct1 <- cor.test(df$wing_mass,  df$mass)
  ct2 <- cor.test(df$wing2_mass, df$mass)
  bind_rows(
    tibble(Metric = "Ratio",  n = nrow(df), r = as.numeric(ct1$estimate), p_value = ct1$p.value),
    tibble(Metric = "Ratio2", n = nrow(df), r = as.numeric(ct2$estimate), p_value = ct2$p.value)
  )
}) %>% list_rbind(names_to = "species_") %>%
  mutate(Study = "Atlantic birds", species = str_replace_all(species_, "_", " "))

write_csv(Ratio_mass_cor, "Derived/Csv/Atlantic_ratio_mass_cor.csv")

# Run models & extract parms ----------------------------------------------
# Data already filtered (unknowns + small groups removed) upstream in Atl_birds_l3.
# Approach-specific covariate rules (Nightjar-aligned):
#   Ratio/Ratio2/Sli_iso: no age/sex by design
#   Resid_ols: age/sex cleaned in first model (Atl_birds_l3); no additional covariates
#   Ryding: include age/sex in combined model (B.Tavg conditional on both)
#   Sli_est: per-group SMA slopes handled upstream in sliR::calc_sli(control = covs)
parms_df <- map(Atl_birds_l4, \(df) {
  sp       <- unique(df$species_)
  covs     <- c(if (sp %in% sig_age_any) "age", if (sp %in% sig_sex_any) "sex")
  covs     <- covs[vapply(covs, \(v) length(unique(na.omit(df[[v]]))) >= 2, logical(1))]
  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  mod_coef_ratio  <- lm(wing_mass    ~ B.Tavg, data = df) %>%
    tidy() %>% mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass   ~ B.Tavg, data = df) %>%
    tidy() %>% mutate(Approach = "Ratio2")
  mod_sli_iso     <- lm(sli_isometry ~ B.Tavg, data = df) %>%
    tidy() %>% mutate(Approach = "Sli_iso")
  mod_resid_ols   <- lm(resid_ols    ~ B.Tavg, data = df) %>%
    tidy() %>% mutate(Approach = "Resid_ols")
  mod_coef_ols    <- lm(as.formula(paste("wing ~ mass + B.Tavg", covs_str)), data = df) %>%
    tidy() %>% mutate(Approach = "Ryding")
  mod_sli_est <- if (any(!is.na(df$sli_estimated))) {
    lm(sli_estimated ~ B.Tavg, data = df) %>% tidy() %>% mutate(Approach = "Sli_est")
  } else tibble()

  bind_rows(mod_coef_ratio, mod_coef_ratio2, mod_coef_ols, mod_resid_ols, mod_sli_est, mod_sli_iso)
}) %>% list_rbind(names_to = "species_") %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

# Rank order summary -------------------------------------------------------
# Groups: Ratio (Ratio+Ratio2), SLI (Sli_est+Sli_iso), OLS (Ryding+Resid_ols)
# Simulation prediction for Bergmann's / Stable:  avg_ratio > avg_sli > avg_ols
# Simulation prediction for Inverse Bergmann's:   avg_ratio < avg_sli < avg_ols
rank_summary <- parms_df %>%
  filter(term == "B.Tavg",
         Approach %in% c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso")) %>%
  mutate(approach_group = case_when(
    Approach %in% c("Ratio", "Ratio2")     ~ "avg_ratio",
    Approach %in% c("Sli_est", "Sli_iso")  ~ "avg_sli",
    Approach %in% c("Ryding", "Resid_ols") ~ "avg_ols"
  )) %>%
  group_by(species_, approach_group) %>%
  summarise(avg_coef = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = approach_group, values_from = avg_coef) %>%
  left_join(Spp_keep2 %>% dplyr::select(species_, Direction), by = "species_") %>%
  mutate(
    rank_consistent = case_when(
      Direction == "Bergmann's"         ~ avg_ratio > avg_sli & avg_sli > avg_ols,
      Direction == "Inverse Bergmann's" ~ avg_ratio < avg_sli & avg_sli < avg_ols,
      TRUE ~ NA
    )
  ) %>%
  arrange(Direction, species_)

rank_summary
rank_summary %>%
  filter(!is.na(rank_consistent)) %>%
  group_by(Direction) %>%
  summarise(n = n(), n_consistent = sum(rank_consistent),
            pct_consistent = round(100 * mean(rank_consistent)))

# Plot slope estimates ----------------------------------------------------
## Slope estimates of temperature's impact on shape

approach_labels <- c(
  "Ratio"    = "Wing / mass",
  "Ratio2"   = "Wing² / mass",
  "Ryding"   = "Mass as covariate",
  "Resid_ols"= "OLS residuals",
  "Sli_est"  = "SLI estimated",
  "Sli_iso"  = "SLI isometry"
)

Num_spp <- length(unique(parms_df$species_))

parms_df_p <- parms_df %>%
  filter(term == "B.Tavg") %>%
  left_join(Spp_metadata2[, c("species_", "Direction", "Sig_trait")]) %>%
  left_join(rank_summary %>% dplyr::select(species_, rank_consistent), by = "species_") %>%
  mutate(Approach = factor(Approach, levels = c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso")),
         species  = str_replace(species_, "_", "\n")) %>%
  filter(std.error < .5)

shape_scale <- c("Both" = 15, "Mass" = 16, "Neither" = 17, "Wing" = 18, "Tarsus" = 19)

direction_order  <- c("Bergmann's", "Inverse Bergmann's", "Mixed - Wingier", "Mixed - Fatter", "Stable")
Direction_effect <- intersect(direction_order, unique(Spp_metadata2$Direction))
Direction_effect <- setNames(Direction_effect, Direction_effect)

plot_shape <- function(df, title = NULL, legend = TRUE, drop_y = FALSE) {
  star_df <- df %>%
    group_by(species) %>%
    summarise(y_star = max(UCI95, na.rm = TRUE),
              rank_consistent = first(rank_consistent),
              .groups = "drop") %>%
    filter(rank_consistent == TRUE)

  p <- df %>%
    mutate(Sig_trait = str_to_sentence(Sig_trait)) %>%
    ggplot(aes(x = species, y = estimate, color = Approach,
               group = interaction(species, Approach))) +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                  alpha = .8, width = 0,
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(shape = Sig_trait),
               size = 2, position = position_dodge(width = 0.75)) +
    geom_text(data = star_df, aes(x = species, y = y_star, label = "*"),
              inherit.aes = FALSE, size = 5, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = NULL, y = expression(beta[T] ~ "on wingyness"),
         title = title) +
    theme(
      axis.text.x = element_text(vjust = .58, angle = 60),
      legend.position = "top"
    ) +
    scale_shape_manual(values = shape_scale, na.value = 16) +
    scale_color_discrete(labels = approach_labels)

  if (drop_y) p <- p + theme(axis.title.y = element_blank())
  if (!legend) p <- p + theme(legend.position = "none")
  return(p)
}

Shape_plots <- imap(Direction_effect, \(direction, name) {
  parms_filt <- parms_df_p %>% filter(Direction == direction)
  drop_y     <- direction %in% c("Inverse Bergmann's", "Mixed - Wingier", "Mixed - Fatter")
  plot_shape(df = parms_filt, title = name, drop_y = drop_y) 
})

# Plot
wrap_plots(Shape_plots, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")

# Save
ggsave("Figures/Atlantic_birds_shape.png", bg = "white", height = 7, width = 11)

# CSV export ---------------------------------------------------------------
write_csv(
  parms_df_p %>%
    mutate(Study = "Atlantic birds", species = str_replace_all(species_, "_", " ")),
  "Derived/Csv/Atlantic_parms.csv"
)

## Interpretation
Spp_metadata2 %>% filter(species_ %in% Spp_keep_vec) %>% 
  mutate(across(.cols = -c(where(is.character), p_mw), round, 2))

# Sicalis flaveola gets 'much wingyer' (i.e. mass decreases and wing increases), so the results are fairly consistent with our simulation
# The species exhibiting inverse Bergmann's rule show the opposite pattern, which makes sense 

# Violin plot
if(FALSE){
  parms_df %>% filter(term == "B.Tavg") %>% 
    ggplot(aes(x = Approach, y = estimate)) +
    geom_violin() + 
    geom_line(aes(group = species_, color = species_)) +
    geom_point(aes(color = species_)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(x = NULL, y = expression(beta[T] ~ "on wingyness"),
         title = paste(Num_spp, "species in analysis")) +
    guides(color = "none") 
}

# Coefficient of variation ------------------------------------------------

calc_lambda <- function(x, y){
  cv_x <- sd({{ x }}) / mean({{ x }})
  cv_y <- sd({{ y }}) / mean({{ y }})
  (cv_y^2) / (cv_x^2)
}

# Very low CVs, suggesting that the variance in wing is very low compared to the variance in mass (relative to the means)
Cvs <- Atlantic_birds5 %>%
  summarize(N = n(), 
            lambda = calc_lambda(x = mass, y = wing),
            .by = species_) %>%
  mutate(lambda = round(lambda, 2)) %>% 
  arrange(lambda) 

# stop()  # exploration checkpoint — commented out to allow full source()

# Interpretation -------------------------------------------------------------

# This doesn't seem to be a great dataset to test this particular question because there are very few species that follow Bergmann's rule and have a significant correlation between wing and mass. Joliceur 1990 argues that the correlation should be greater than 0.6, which is extremely high. Furthermore, the relationships between temperature and size are pretty small, so temperature may have little effect on size and shape.
#NOTE: I tried with tarsus as well and got similar results

# Next steps --------------------------------------------------------------

# Could try this looking at how shape changes with the amount of forest or the habitat where the bird was captured
# Try with year instead of average temperature. Climate change / deforestation (also increasing through time) may be stronger influences than temperature changes by geography
# Could try with some sort of bill measurement

# EXTRAS ------------------------------------------------------------------

Spp_metadata2 %>% 
  #filter(cor_mw > 0.3) %>%
  ggplot(aes(x = b_ols, y = b_sma, color = cor_mw, shape = Keep)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

# Test effect of approach -------------------------------------------------
# stop()  # exploration checkpoint
# Does CV have an impact on the estimate?
if(FALSE){
  parms_df %>% filter(term == "B.Tavg") %>% 
    full_join(Cvs) %>% 
    ggplot(aes(x = lambda, y = estimate, color = Approach)) + 
    geom_point() + 
    geom_smooth(se = FALSE) 
}

parms_df %>% filter(term == "B.Tavg") %>% 
  full_join(Cvs) %>% 
  ggplot(aes(x = Approach, y = estimate, color = species_)) + 
  geom_point() + 
  geom_line(aes(group = species_)) + 
  guides(color = "none")

## GPT
df <- parms_df %>% 
  filter(term == "B.Tavg") %>% 
  full_join(Cvs)

# reorder approaches by mean estimate
df <- df %>% 
  group_by(Approach) %>% 
  mutate(mean_est = mean(estimate, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Approach = reorder(Approach, mean_est))

# overall mean trend for black line
overall <- df %>% 
  group_by(Approach) %>% 
  summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = "drop")

ggplot(df, aes(x = Approach, y = estimate, color = species_)) + 
  geom_point() +
  geom_line(aes(group = species_), alpha = 0.5) +
  guides(color = "none") +
  # black trend line of overall mean
  geom_line(data = overall, aes(x = Approach, y = mean_est, group = 1), 
            color = "black", size = 1) +
  theme_minimal(base_size = 14)

# Test assumptions --------------------------------------------------------
# (leftover nightjar code — nj_df_l does not exist in this script)
# ols_mod_l <- map(nj_df_l, \(df){ lm(log_wing ~ log_mass + B.Tavg, data = df) })
# sma_mod_l <- map(nj_df_l, \(df){ sma(log_wing ~ log_mass, data = df, method = "SMA") })
