## Temporal shape-shifting in North American migratory passerines ##
## Data: Weeks et al. (2020) Ecology Letters — 70,716 capture records, 53 species, 1979-2016

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(smatr)
library(cowplot)
library(broom)
library(janitor)
library(ggpubr)
ggplot2::theme_set(theme_cowplot())

source("Scripts/Key_allometry_fns.R")

Weeks_path <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Weeks_etal_2020_Data.csv"
Weeks20    <- read_csv(Weeks_path, skip = 2)

# Run parameters ----------------------------------------------------------
min_n_obs       <- 150    # minimum observations per species
year_cutoff     <- 1978   # exclude records from this year and before
p_bergmann      <- 0.05   # threshold for temporal trend classification (mass, wing)
p_age_sex       <- 0.10   # threshold for age / sex covariate inclusion
pos_allom       <- TRUE   # retain only species with b_ols > cor_b_min & p_mw < cor_p_max
filter_temporal <- FALSE  # TRUE: retain only species with a significant temporal mass trend
control_age_sex <- TRUE   # TRUE: include age / sex covariates where they significantly affect morphology
min_n_age_group <- 100    # min individuals per age group (AHY/HY) to include age control
min_n_sex_group <- 100    # min individuals per sex group (m/f) to include sex control
cor_b_min       <- 0.1    # min b_ols (Mass ~ Wing) within group for allometric correlation filter
cor_p_max       <- 0.05   # max p-value for Mass ~ Wing within group

# Format ------------------------------------------------------------------
# NOTE on sign convention: positive Year coefficient = trait increases over time.
# Weeks 2020 finding: Mass decreases (negative Year coef) and Wing increases (positive),
# i.e. birds are becoming relatively longer-winged. We test whether this holds after
# accounting for latent body size.

Weeks_df <- Weeks20 %>%
  rename(species_ = Taxon) %>%
  mutate(year       = as.integer(format(mdy(Date), "%Y")),
         log_mass   = log(Mass),
         log_wing   = log(Wing),
         log_tarsus = log(Tarsus),
         wing_mass  = log_wing - log_mass,
         wing2_mass = 2*log_wing - log_mass) %>%
  dplyr::select(species_, Mass, Wing, Tarsus, year, Sex, Age,
                log_mass, log_wing, log_tarsus, wing_mass, wing2_mass) %>%
  filter(!is.na(Mass) & !is.na(Wing) & !is.na(Tarsus)) %>%
  filter(year > year_cutoff)

# How many species and complete-case observations?
Weeks_df %>% summarize(n = n(), .by = species_) %>% arrange(desc(n))

# Species-level filtering -------------------------------------------------
Spp_summary <- Weeks_df %>%
  summarize(num_obs    = n(),
            min_year   = min(year),
            max_year   = max(year),
            year_range = max_year - min_year,
            .by = species_)

Spp_include <- Spp_summary %>%
  filter(num_obs >= min_n_obs) %>%
  arrange(desc(num_obs))
Spp_include

# Create per-species list
Weeks_l <- Weeks_df %>%
  right_join(Spp_include %>% dplyr::select(species_), by = "species_") %>%
  group_split(species_)
names(Weeks_l) <- map_chr(Weeks_l, ~unique(pull(.x, species_)))
length(Weeks_l)

# Age & Sex ---------------------------------------------------------------
Weeks_df %>% tabyl(Age)
Weeks_df %>% tabyl(Sex)

# Test whether Age / Sex significantly affect each morphometric, controlling for year.
# Age levels: AHY (after-hatch-year = adult) and HY (hatch-year = juvenile).
# Sex levels: m, f.
test_group_effect <- function(df_list, dv, iv, temp_name = "year",
                              min_n_per_group = 0) {
  map(df_list, \(df) {
    df_filt <- df %>% filter(!is.na(.data[[iv]]) & .data[[iv]] != "Unknown")
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
    fmla <- as.formula(paste(dv, "~", temp_name, "+", iv))
    tidy(lm(fmla, data = df_filt)) %>%
      filter(str_starts(term, iv)) %>%
      mutate(dv = dv) %>%
      bind_cols(n_counts)
  }) %>%
    list_rbind(names_to = "species_")
}

Age_tbl <- bind_rows(
  test_group_effect(Weeks_l, dv = "Mass",   iv = "Age", min_n_per_group = min_n_age_group),
  test_group_effect(Weeks_l, dv = "Wing",   iv = "Age", min_n_per_group = min_n_age_group),
  test_group_effect(Weeks_l, dv = "Tarsus", iv = "Age", min_n_per_group = min_n_age_group)
) %>% mutate(sig = p.value < p_age_sex)

Sex_tbl <- bind_rows(
  test_group_effect(Weeks_l, dv = "Mass",   iv = "Sex", min_n_per_group = min_n_sex_group),
  test_group_effect(Weeks_l, dv = "Wing",   iv = "Sex", min_n_per_group = min_n_sex_group),
  test_group_effect(Weeks_l, dv = "Tarsus", iv = "Sex", min_n_per_group = min_n_sex_group)
) %>% mutate(sig = p.value < p_age_sex)

Age_tbl %>% filter(sig) %>% dplyr::select(species_, dv, estimate, p.value)
Sex_tbl %>% filter(sig) %>% dplyr::select(species_, dv, estimate, p.value)

sig_age_mass   <- Age_tbl %>% filter(dv == "Mass",   sig) %>% pull(species_)
sig_sex_mass   <- Sex_tbl %>% filter(dv == "Mass",   sig) %>% pull(species_)
sig_age_wing   <- Age_tbl %>% filter(dv == "Wing",   sig) %>% pull(species_)
sig_sex_wing   <- Sex_tbl %>% filter(dv == "Wing",   sig) %>% pull(species_)
sig_age_tarsus <- Age_tbl %>% filter(dv == "Tarsus", sig) %>% pull(species_)
sig_sex_tarsus <- Sex_tbl %>% filter(dv == "Tarsus", sig) %>% pull(species_)

# Union: if age/sex affects any morphometric, include in all downstream models
sig_age_any <- unique(c(sig_age_mass, sig_age_wing, sig_age_tarsus))
sig_sex_any <- unique(c(sig_sex_mass, sig_sex_wing, sig_sex_tarsus))

if (!control_age_sex) {
  sig_age_any <- character(0)
  sig_sex_any <- character(0)
}

# Temporal mass / wing trend (Bergmann's analog) --------------------------
# Per species: does Mass or Wing change significantly over time?
# Weeks 2020: mass shrinks (negative coef), wing grows (positive coef).

temporal_model <- function(df, dv, sig_age_spp, sig_sex_spp) {
  sp   <- unique(df$species_)
  df_  <- df
  covs <- character(0)
  if (sp %in% sig_age_spp) {
    df_  <- df_ %>% filter(!is.na(Age))
    covs <- c(covs, "Age")
  }
  if (sp %in% sig_sex_spp) {
    df_  <- df_ %>% filter(!is.na(Sex))
    covs <- c(covs, "Sex")
  }
  rhs <- paste(c("year", covs), collapse = " + ")
  tidy(lm(as.formula(paste(dv, "~", rhs)), data = df_)) %>%
    filter(term == "year") %>%
    mutate(dv = dv)
}

mass_mod <- map(Weeks_l, temporal_model, dv = "Mass",
                sig_age_spp = sig_age_mass, sig_sex_spp = sig_sex_mass) %>%
  list_rbind(names_to = "species_")

wing_mod <- map(Weeks_l, temporal_model, dv = "Wing",
                sig_age_spp = sig_age_wing, sig_sex_spp = sig_sex_wing) %>%
  list_rbind(names_to = "species_")

Temporal_mods <- bind_rows(mass_mod, wing_mod)

# Classify direction and significance of temporal trends
Spp_keep <- Temporal_mods %>%
  mutate(Increasing = ifelse(estimate > 0, "Y", "N")) %>%
  distinct(species_, Increasing, dv, p.value) %>%
  arrange(species_)

Spp_keep2 <- Spp_keep %>%
  mutate(sig = p.value < p_bergmann) %>%
  group_by(species_) %>%
  summarise(
    n_sig     = sum(sig),
    mass_dir  = Increasing[dv == "Mass"],   # "Y" = increasing over time, "N" = decreasing
    wing_dir  = Increasing[dv == "Wing"],
    mass_sig  = sig[dv == "Mass"],
    wing_sig  = sig[dv == "Wing"],
    Sig_trait = ifelse(n_sig == 0, "Neither",
                       ifelse(n_sig == 2, "both", dv[sig])),
    .groups   = "drop"
  ) %>%
  mutate(
    sole_dir  = case_when(mass_sig & !wing_sig ~ mass_dir,
                          !mass_sig & wing_sig ~ wing_dir),
    Direction = case_when(
      n_sig == 0                                          ~ "Stable",
      n_sig == 1 & sole_dir == "N"                       ~ "Bergmann's",
      n_sig == 1 & sole_dir == "Y"                       ~ "Inverse Bergmann's",
      n_sig == 2 & mass_dir == "N" & wing_dir == "N"     ~ "Bergmann's",
      n_sig == 2 & mass_dir == "Y" & wing_dir == "Y"     ~ "Inverse Bergmann's",
      n_sig == 2 & mass_dir == "N" & wing_dir == "Y"     ~ "Mixed - Wingier",
      n_sig == 2 & mass_dir == "Y" & wing_dir == "N"     ~ "Mixed - Fatter",
      TRUE ~ "Check"
    )
  )

Spp_keep2 %>% tabyl(Direction)

# Visualize mass and wing temporal trends per species
Temporal_mods %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error) %>%
  ggplot(aes(x = reorder(species_, estimate), y = estimate,
             color = dv, ymin = LCI95, ymax = UCI95)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point() +
  labs(x = NULL, y = "Year coefficient (raw scale)", color = NULL) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

if (filter_temporal) Spp_keep2 <- Spp_keep2 %>% filter(Direction != "Stable")

Weeks_df2 <- Weeks_df %>% filter(species_ %in% Spp_keep2$species_)

# Allometric correlation --------------------------------------------------
# SMA is appropriate only when mass and wing are meaningfully correlated.
# Filter to species with positive allometric slope and significant Wing-Mass correlation.

Slopes_tbl <- Weeks_df2 %>%
  summarize(b_sma = sd(log_wing, na.rm = TRUE) / sd(log_mass, na.rm = TRUE),
            b_ols = cor(Wing, Mass, use = "complete.obs") * b_sma,
            .by = species_)

Spp_metadata <- Spp_keep2 %>%
  left_join(Spp_include, by = "species_") %>%
  left_join(Slopes_tbl,  by = "species_")
Spp_metadata

Weeks_l2 <- Weeks_df2 %>% group_split(species_)
names(Weeks_l2) <- map_chr(Weeks_l2, ~unique(pull(.x, species_)))

Cors_tbl <- map(Weeks_l2, \(df) {
  ct <- cor.test(df$Mass, df$Wing, use = "complete.obs")
  tibble(cor_mw = as.numeric(ct$estimate), p_mw = ct$p.value)
}) %>% list_rbind(names_to = "species_")

Spp_metadata2 <- Spp_metadata %>% left_join(Cors_tbl) %>% 
  mutate(Keep = ifelse(b_ols > cor_b_min & p_mw < cor_p_max, "Include", "Exclude"))

Spp_keep_vec <- Spp_metadata2 %>%
  filter(Keep == "Include") %>%
  pull(species_)
Spp_keep_vec

# Weeks_l2 is NOT filtered to Spp_keep_vec here; sli_est is set NA for non-passing species in Weeks_l3
if (pos_allom) message(length(Spp_keep_vec), " / ", length(Weeks_l2), " species pass allometric filter (sli_est will be NA for the rest)")

# Visualize allometry across retained species
Weeks_df2 %>%
  filter(species_ %in% names(Weeks_l2)) %>%
  ggplot(aes(x = log_mass, y = log_wing, color = species_)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  ggpmisc::stat_ma_line(method = "SMA", se = FALSE) +
  guides(color = "none") +
  labs(title = paste(length(Weeks_l2), "species after allometric filter"))

# Per-group allometric correlation (Mass ~ Wing within each Age × Sex combination) -----
if (control_age_sex) {
  group_cor_wing <- imap(Weeks_l2, \(df, sp) {
    covs <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
    if (!length(covs)) return(NULL)
    build_group_cor_tbl(df, Append = Wing, Mass = Mass, control = covs) %>%
      mutate(species_ = sp, .before = 1)
  }) %>% list_rbind() %>%
    mutate(pass = b_ols > cor_b_min & p_mw <= cor_p_max)
  print(group_cor_wing)
}

# Shape metrics -----------------------------------------------------------
# Per species: filter to groups meeting min_n thresholds; include Age/Sex in
# covariate list only if ≥ 2 valid groups remain after filtering.
Weeks_l3 <- imap(Weeks_l2, \(df, sp) {
  covs <- character(0)

  if (sp %in% sig_age_any) {
    valid_age <- df %>%
      filter(!is.na(Age)) %>%
      count(Age) %>% filter(n >= min_n_age_group) %>% pull(Age)
    if (length(valid_age) >= 1) df <- df %>% filter(Age %in% valid_age)
    if (length(valid_age) >= 2) {
      passing_age <- build_group_cor_tbl(df, Append = Wing, Mass = Mass, control = "Age") %>%
        ungroup() %>% filter(b_ols > cor_b_min & p_mw <= cor_p_max) %>% pull(Age)
      if (length(passing_age) >= 1) df <- df %>% filter(Age %in% passing_age)
      if (length(passing_age) >= 2) covs <- c(covs, "Age")
    }
  }
  if (sp %in% sig_sex_any) {
    valid_sex <- df %>%
      filter(!is.na(Sex)) %>%
      count(Sex) %>% filter(n >= min_n_sex_group) %>% pull(Sex)
    if (length(valid_sex) >= 1) df <- df %>% filter(Sex %in% valid_sex)
    if (length(valid_sex) >= 2) {
      passing_sex <- build_group_cor_tbl(df, Append = Wing, Mass = Mass, control = "Sex") %>%
        ungroup() %>% filter(b_ols > cor_b_min & p_mw <= cor_p_max) %>% pull(Sex)
      if (length(passing_sex) >= 1) df <- df %>% filter(Sex %in% passing_sex)
      if (length(passing_sex) >= 2) covs <- c(covs, "Sex")
    }
  }

  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  ols_wing   <- lm(as.formula(paste("log_wing   ~ log_mass", covs_str)), data = df)
  ols_tarsus <- lm(as.formula(paste("log_tarsus ~ log_mass", covs_str)), data = df)
  sma_wing   <- sma(log_wing  ~ log_mass, data = df, method = "SMA")
  sma_tarsus <- sma(log_tarsus ~ log_mass, data = df, method = "SMA")

  df_res <- df %>%
    mutate(resid_ols        = log_wing   - predict(ols_wing,   newdata = df),
           resid_sma        = residuals(sma_wing),
           resid_ols_tarsus = log_tarsus - predict(ols_tarsus, newdata = df),
           resid_sma_tarsus = residuals(sma_tarsus)) %>%
    calc_sli(b_sli = 0.33, Append = Wing,   rename_col = "sli_isometry") %>%
    calc_sli(b_sli = 0.33, Append = Tarsus, rename_col = "sli_tarsus_iso")

  if (sp %in% Spp_keep_vec) {
    df_res %>%
      calc_sli(Append = Wing,
               control = if (length(covs)) covs else NULL,
               b_sli   = coef(sma_wing)["slope"],   rename_col = "sli_estimated") %>%
      calc_sli(Append = Tarsus,
               control = if (length(covs)) covs else NULL,
               b_sli   = coef(sma_tarsus)["slope"], rename_col = "sli_tarsus_est")
  } else {
    df_res %>% mutate(sli_estimated = NA_real_, sli_tarsus_est = NA_real_)
  }
})

# Per-group SMA slope table (for inspection — wing)
if (control_age_sex) {
  wing_slopes_tbl <- imap(Weeks_l2, \(df, sp) {
    covs <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
    if (!length(covs)) return(NULL)
    build_sli_slopes_tbl(df, Append = Wing, Mass = Mass, control = covs) %>%
      mutate(species_ = sp, .before = 1)
  }) %>% list_rbind()
  print(wing_slopes_tbl)
}

# Z-score per species (required for standardised comparisons)
Weeks_l4 <- map(Weeks_l3, \(df){
  df %>% mutate(across(where(is.numeric), scale))
})

# Inspect correlation of shape metrics with body size (mass)
map(Weeks_l4, \(df){
  df %>% summarize(wm_m    = cor(wing_mass, Mass),
                   w2m_m   = cor(wing2_mass, Mass),
                   resid_m = cor(resid_sma, Mass))
})

# Run models & extract parms ----------------------------------------------
# Data already filtered (small groups removed) upstream in Weeks_l3.
# Approach-specific covariate rules (Nightjar-aligned):
#   Ratio/Ratio2/Sli_iso: no Age/Sex by design
#   Resid_ols: Age/Sex cleaned in first model (Weeks_l3); no additional covariates
#   Ryding: include Age/Sex in combined model (year conditional on both)
#   Sli_est: per-group SMA slopes handled upstream in calc_sli(control = covs)
parms_df <- map(Weeks_l4, \(df) {
  sp       <- unique(df$species_)
  covs     <- c(if (sp %in% sig_age_any) "Age", if (sp %in% sig_sex_any) "Sex")
  covs     <- covs[vapply(covs, \(v) length(unique(na.omit(df[[v]]))) >= 2, logical(1))]
  covs_str <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  mod_coef_ratio  <- lm(wing_mass    ~ year, data = df) %>%
    tidy() %>% mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass   ~ year, data = df) %>%
    tidy() %>% mutate(Approach = "Ratio2")
  mod_sli_iso     <- lm(sli_isometry ~ year, data = df) %>%
    tidy() %>% mutate(Approach = "Sli_iso")
  mod_resid_ols   <- lm(resid_ols    ~ year, data = df) %>%
    tidy() %>% mutate(Approach = "Resid_ols")
  mod_coef_ols    <- lm(as.formula(paste("Wing ~ Mass + year", covs_str)), data = df) %>%
    tidy() %>% mutate(Approach = "Ryding")
  mod_sli_est <- if (any(!is.na(df$sli_estimated))) {
    lm(sli_estimated ~ year, data = df) %>% tidy() %>% mutate(Approach = "Sli_est")
  } else tibble()

  bind_rows(mod_coef_ratio, mod_coef_ratio2, mod_coef_ols, mod_resid_ols, mod_sli_est, mod_sli_iso)
}) %>% list_rbind(names_to = "species_") %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

# Rank order summary -------------------------------------------------------
# Groups: Ratio (Ratio+Ratio2), SLI (Sli_est+Sli_iso), OLS (Ryding+Resid_ols)
# Simulation prediction for Bergmann's species:  avg_ratio > avg_sli > avg_ols
# Simulation prediction for Inverse Bergmann's:  avg_ratio < avg_sli < avg_ols
rank_summary <- parms_df %>%
  filter(term == "year",
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
approach_labels <- c(
  "Ratio"      = "Wing / Mass",
  "Ratio2"     = "Wing² / Mass",
  "Ryding"     = "Mass as covariate",
  "Resid_ols"  = "OLS residuals",
  "Sli_est"    = "SLI estimated",
  "Sli_iso"    = "SLI isometry"
)

Num_spp <- length(unique(parms_df$species_))

parms_df_p <- parms_df %>%
  filter(term == "year") %>%
  left_join(Spp_metadata2[, c("species_", "Direction", "Sig_trait")], by = "species_") %>%
  left_join(rank_summary %>% dplyr::select(species_, rank_consistent), by = "species_") %>%
  mutate(Approach = factor(Approach, levels = c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso")),
         species  = str_replace(species_, "_", " ")) %>%
  filter(std.error < .5)

direction_order  <- c("Bergmann's", "Inverse Bergmann's", "Mixed - Wingier", "Mixed - Fatter", "Stable")
Direction_effect <- intersect(direction_order, unique(parms_df_p$Direction))
Direction_effect <- setNames(Direction_effect, Direction_effect)

shape_scale <- c("Both" = 15, "Mass" = 16, "Neither" = 17, "Wing" = 18, "Tarsus" = 19)

plot_shape <- function(df, title = NULL, legend = TRUE, drop_y = FALSE) {
  star_df <- df %>%
    group_by(species) %>%
    summarise(y_star = max(UCI95, na.rm = TRUE),
              rank_consistent = first(rank_consistent),
              .groups = "drop") %>%
    filter(rank_consistent == TRUE)

  p <- df %>%
    mutate(Sig_trait = str_to_sentence(Sig_trait)) %>%
    ggplot(aes(x = reorder(species, estimate), y = estimate, color = Approach,
               group = interaction(species, Approach))) +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                  alpha = .8, width = 0,
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(shape = Sig_trait), size = 2, position = position_dodge(width = 0.75)) +
    geom_text(data = star_df, aes(x = species, y = y_star, label = "*"),
              inherit.aes = FALSE, size = 5, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = NULL, y = expression(beta[Year] ~ "on winginess"),
         title = title) +
    theme(axis.text.x = element_text(hjust = 1, angle = 60),
          legend.position = "top") +
    scale_shape_manual(values = shape_scale) +
    scale_color_discrete(labels = approach_labels)
  if (drop_y) p <- p + theme(axis.title.y = element_blank())
  if (!legend) p <- p + theme(legend.position = "none")
  return(p)
}

Shape_plots <- imap(Direction_effect, \(direction, name) {
  parms_filt <- parms_df_p %>% filter(Direction == direction)
  plot_shape(df = parms_filt, title = name)
})

com.leg <- get_legend(plot_shape(df = parms_df_p))
ggarrange(plotlist = Shape_plots, common.legend = TRUE, legend.grob = com.leg, labels = "auto")
ggsave("Figures/Weeks_shape.png", bg = "white", height = 7, width = max(9, Num_spp * 0.4))

# CSV export ---------------------------------------------------------------
write_csv(
  parms_df_p %>% mutate(Study = "Weeks (2020)"),
  "Derived/Csv/Weeks_parms.csv"
)

# EXTRAS ------------------------------------------------------------------

Spp_metadata2 %>% 
  #filter(cor_mw > 0.3) %>%
  ggplot(aes(x = b_ols, y = b_sma, color = cor_mw, shape = Keep)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

# Species that are 'mixed - wingier' but overlap 0 
Spp_metadata2 %>% 
  filter(species_ %in% c("Setophaga_palmarum", "Oreothlypis_ruficapilla", "Oporornis_agilis"))
