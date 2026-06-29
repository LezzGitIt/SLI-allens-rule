## Examining changes in wing shape in the Brazilian Atlantic Forest ##

# Libraries & data -------------------------------------------------
library(tidyverse)
library(broom)
library(janitor)
library(smatr)
library(cowplot)
library(ggpubr)
library(naniar)
library(lavaan)
source("Scripts/Key_causal_fns.R")
ggplot2::theme_set(theme_cowplot())

Atlantic_birds <- read_csv("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Ecology/Atlantic_bird_traits/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv")

# Run parameters ----------------------------------------------------------
min_n_obs       <- 150    # minimum observations per species to be included
year_cutoff     <- 1990   # exclude records from this year and before
p_bergmann      <- 0.05   # p-value threshold: Bergmann / Allen response classification
p_age_sex       <- 0.10   # p-value threshold: age / sex covariate inclusion
pos_allom       <- TRUE   # retain only species with positive allometric slope (b_ols > 0.1, p < .06)
filter_bergmann <- FALSE  # TRUE: retain only species with a significant Bergmann / Allen response
control_age_sex <- TRUE   # TRUE: include age / sex as covariates where they significantly affect body size

# SLI function ------------------------------------------------------------
# Calculate b_sli using logged morphometrics
calc_sli <- function(df, Append = Append, b_sli = 0.33, rename_col = FALSE){
  # L0 is the average mass, essentially allowing for comparison of wing lengths for a given mass
  L0 <- mean(df$mass)
  df_sli <- df %>% mutate(sli = {{ Append }} * (L0 / mass)^b_sli) %>%
    arrange(desc(sli))
  if(rename_col != FALSE){df_sli <- df_sli %>% rename( {{ rename_col }} := sli)}
  return(df_sli)
}

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
         wing_mass = wing / mass, 
         wing2_mass = wing^2 / mass)

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
# Filters Unknown/NA for iv, includes per-group sample counts in the output.
test_group_effect <- function(df_list, dv, iv) {
  map(df_list, \(df) {
    df_filt <- df %>% filter(.data[[iv]] != "Unknown" & !is.na(.data[[iv]]))
    if (nrow(df_filt) < 10 || length(unique(df_filt[[iv]])) < 2) return(NULL)

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
Age_tbl <- bind_rows(
  test_group_effect(Atl_birds_l, dv = "mass",   iv = "age"),
  test_group_effect(Atl_birds_l, dv = "wing",   iv = "age"),
  test_group_effect(Atl_birds_l, dv = "tarsus", iv = "age")
) %>% mutate(sig = p.value < p_age_sex) %>%
  filter(n_juvenile > 20)

Sex_tbl <- bind_rows(
  test_group_effect(Atl_birds_l, dv = "mass",   iv = "sex"),
  test_group_effect(Atl_birds_l, dv = "wing",   iv = "sex"),
  test_group_effect(Atl_birds_l, dv = "tarsus", iv = "sex")
) %>% mutate(sig = p.value < p_age_sex) %>%
  filter(n_female > 20)

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

# Union across all DVs: if age/sex affects any indicator, include it in SEM + shape models
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

# Keep species with significant response ----------------------------------
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
    n_sig = sum(sig),
    Berg_vals = paste0(Bergs[sig], collapse = ""),
    Sig_trait = ifelse(n_sig == 0, "Neither", ifelse(n_sig == 2, "both", dv[sig])),
    .groups = "drop"
  ) %>% 
  mutate(
    Direction = case_when(
      n_sig == 0 ~ "Neither",
      n_sig == 1 & Berg_vals == "Y" ~ "Bergmann's",
      n_sig == 1 & Berg_vals == "N" ~ "Inverse",
      n_sig == 2 & Berg_vals %in% c("YY") ~ "Bergmann's",
      n_sig == 2 & Berg_vals %in% c("NN") ~ "Inverse",
      n_sig == 2 & Berg_vals %in% c("YN", "NY") ~ "Mixed",
      TRUE ~ "Check"
    )
  )
Spp_keep2 %>% tabyl(Direction)

# If filter_bergmann = TRUE, drop species with no significant Bergmann / Allen response
if (filter_bergmann) Spp_keep2 <- Spp_keep2 %>% filter(Direction != "Neither")

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

# Extract the correlations
Cors_tbl <- map(Atl_birds_l2, \(df){
  # In a variable w/ single predictor this is the correlation
  tidy(lm(mass ~ wing, data = df)) 
}) %>% list_rbind(names_to = "species_") %>% 
  filter(term == "wing") %>% 
  dplyr::select(species_, estimate, p.value) %>% 
  rename(cor_mw = estimate, p_mw = p.value)

Spp_metadata2 <- Spp_metadata %>% left_join(Cors_tbl)

# 8 species have significant correlation between mass and wing!!
Spp_keep_vec <- Spp_metadata2 %>% 
  filter(b_ols > 0.1 & p_mw < .06) %>%
  pull(species_)
Spp_keep_vec

# Keep just species with significant correlation
if(pos_allom){
  Atl_birds_l2 <- Atl_birds_l2[Spp_keep_vec]
}

# Create shape metrics ----------------------------------------------------
Atl_birds_l3 <- map(Atl_birds_l2, \(df){
  ols_wing     <- lm(log_wing ~ log_mass, data = df)
  sma_wing     <- sma(log_wing ~ log_mass, data = df, method = "SMA")
  ols_tarsus   <- lm(log_tarsus ~ log_mass, data = df)
  sma_tarsus   <- sma(log_tarsus ~ log_mass, data = df, method = "SMA")
  b_sma_wing   <- coef(sma_wing)['slope']
  b_sma_tarsus <- coef(sma_tarsus)['slope']
  df %>%
    mutate(resid_ols        = residuals(ols_wing),
           resid_sma        = residuals(sma_wing),
           resid_ols_tarsus = residuals(ols_tarsus),
           resid_sma_tarsus = residuals(sma_tarsus)) %>%
    calc_sli(b_sli = 0.33,          Append = wing,   rename_col = "sli_isometry") %>%
    calc_sli(b_sli = b_sma_wing,    Append = wing,   rename_col = "sli_estimated") %>%
    calc_sli(b_sli = 0.33,          Append = tarsus, rename_col = "sli_tarsus_iso") %>%
    calc_sli(b_sli = b_sma_tarsus,  Append = tarsus, rename_col = "sli_tarsus_est")
})

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

# Run models & extract parms ----------------------------------------------
# For each species, age/sex dummies are added as covariates when they significantly affect any of mass, wing, or tarsus (sig_age_any / sig_sex_any).
# Unknown/NA observations are filtered only for species where the covariate matters.
parms_df <- map(Atl_birds_l4, \(df) {
  sp   <- unique(df$species_)
  df_  <- df
  covs <- character(0)
  if (sp %in% sig_age_any) {
    df_  <- df_ %>% filter(age != "Unknown" & !is.na(age))
    covs <- c(covs, "age")
  }
  if (sp %in% sig_sex_any) {
    df_  <- df_ %>% filter(sex != "Unknown" & !is.na(sex))
    covs <- c(covs, "sex")
  }
  extra <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  mod_resid_ols   <- lm(as.formula(paste("resid_ols ~ B.Tavg",        extra)), data = df_) %>%
    tidy() %>% mutate(Approach = "Resid_ols")
  mod_coef_ols    <- lm(as.formula(paste("wing ~ mass + B.Tavg",      extra)), data = df_) %>%
    tidy() %>% mutate(Approach = "Ryding")
  mod_coef_ratio  <- lm(as.formula(paste("wing_mass ~ B.Tavg",        extra)), data = df_) %>%
    tidy() %>% mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(as.formula(paste("wing2_mass ~ B.Tavg",       extra)), data = df_) %>%
    tidy() %>% mutate(Approach = "Ratio2")
  mod_sli_iso     <- lm(as.formula(paste("sli_isometry ~ B.Tavg",     extra)), data = df_) %>%
    tidy() %>% mutate(Approach = "Sli_iso")
  mod_sli_est     <- lm(as.formula(paste("sli_estimated ~ B.Tavg",    extra)), data = df_) %>%
    tidy() %>% mutate(Approach = "Sli_est")
  bind_rows(mod_sli_iso, mod_sli_est, mod_resid_ols, mod_coef_ols, mod_coef_ratio, mod_coef_ratio2)
}) %>% list_rbind(names_to = "species_") %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

# Add SEM: Atl_birds_l4 is already z-scored per species.
# When age/sex affects any indicator, the column name is passed to fit_lavaan_sem() via
# size_covs; the function handles numeric encoding internally (lavaan requires numeric).
sem_parms <- map(Atl_birds_l4, \(df) {
  sp        <- unique(df$species_)
  df_       <- df
  size_covs <- character(0)
  if (sp %in% sig_age_any) {
    df_       <- df_ %>% filter(age != "Unknown" & !is.na(age))
    size_covs <- c(size_covs, "age")
  }
  if (sp %in% sig_sex_any) {
    df_       <- df_ %>% filter(sex != "Unknown" & !is.na(sex))
    size_covs <- c(size_covs, "sex")
  }

  res <- fit_lavaan_sem(
    df_,
    mass_name    = "log_mass",
    append_names = c("log_wing", "log_tarsus"),
    temp_name    = "B.Tavg",
    labels       = c("wing", "tarsus"),
    size_covs    = size_covs
  )
  tibble(
    term      = "B.Tavg",
    estimate  = c(res$coef_sem_wing,       res$coef_sem_tarsus),
    std.error = c(res$se_sem_wing,         res$se_sem_tarsus),
    lambda    = c(res$lambda_sem_wing,     res$lambda_sem_tarsus),
    se_lambda = c(res$se_lambda_sem_wing,  res$se_lambda_sem_tarsus),
    p.value   = NA_real_,
    Approach  = c("SEM_wing", "SEM_tarsus")
  )
}) %>%
  list_rbind(names_to = "species_") %>%
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

parms_df <- bind_rows(parms_df, sem_parms)

# SEM diagnostics: factor loadings  ---------------------------------------
# lambda = loading of each appendage on latent Size.
# It captures how much of that appendage's variance is shared with mass.
#
# Boundary cases:
#   low lambda  → appendage barely tracks body size; little collider bias to
#                 correct, but the latent factor is poorly constrained (high SE).
#   high lambda → appendage is dominated by size variation; the direct Temp
#                 effect has little independent variance left (also high SE).
#   optimal     → intermediate lambda gives best precision.
#
# If the SEM guard (|coef| > 1.5) discards a species, lambda is NA.

sem_diag <- sem_parms %>%
  filter(!is.na(lambda)) %>%
  left_join(Spp_keep2 %>% dplyr::select(species_, Direction), by = "species_") %>%
  mutate(
    appendage = if_else(Approach == "SEM_wing", "Wing", "Tarsus"),
    species   = str_replace(species_, "_", " "),
    ci_width  = UCI95 - LCI95
  )

# Lambda vs SE of temperature coefficient
p_lambda_se <- ggplot(sem_diag, aes(x = lambda, y = std.error, color = appendage)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  geom_smooth(se = FALSE, method = "loess", span = 1.5, linewidth = 0.8) +
  scale_color_manual(values = c("Wing" = "steelblue", "Tarsus" = "coral3")) +
  annotate("text", x = -Inf, y = 0.52, hjust = -0.1, size = 3,
           label = "SE = 0.5", color = "grey40") +
  labs(x = expression(lambda ~ "(factor loading on latent Size)"),
       y = "SE of temperature coefficient",
       color = NULL)

# Lambda vs coefficient estimate (with 95% CI error bars)
p_lambda_coef <- ggplot(sem_diag, aes(x = lambda, y = estimate, color = appendage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = LCI95, ymax = UCI95), width = 0, alpha = 0.35) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Wing" = "steelblue", "Tarsus" = "coral3")) +
  labs(x = expression(lambda ~ "(factor loading on latent Size)"),
       y = "Temperature coefficient",
       color = NULL)

# Lambda and its own SE: species where the loading is imprecise
p_lambda_uncertainty <- ggplot(sem_diag, aes(x = lambda, y = se_lambda, color = appendage)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Wing" = "steelblue", "Tarsus" = "coral3")) +
  labs(x = expression(lambda ~ "(factor loading on latent Size)"),
       y = expression("SE of" ~ lambda),
       color = NULL)

plot_grid(p_lambda_se, p_lambda_coef, p_lambda_uncertainty, nrow = 1)

# Strong relationship between the uncertainty in factor loadings and the uncertainty in the temperature coefficient 
sem_diag %>% filter(se_lambda < 10) %>%
  ggplot(aes(x = se_lambda, y = std.error, color = appendage)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Wing" = "steelblue", "Tarsus" = "coral3")) +
  labs(x = expression("SE of" ~ lambda),
       y = "SE of temperature coefficient",
       color = NULL)

# Not a strong relationship with sample size and SE of lambda, although all the ones with very high lambda SEs have low sample size
sem_diag %>% #filter(se_lambda < 5) %>%
  left_join(Spp_metadata2[,c("species_", "num_obs")]) %>%
  ggplot(aes(x = num_obs, y = se_lambda, color = appendage)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Wing" = "steelblue", "Tarsus" = "coral3")) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = "Sample size",
       y = expression("SE of" ~ lambda),
       color = NULL)

# Tabulate: which species have high-SE SEM estimates?
sem_diag %>%
  dplyr::select(species_, appendage, lambda, se_lambda, estimate, std.error, LCI95, UCI95, Direction) %>%
  arrange(appendage, desc(se_lambda))

# Plot slope estimates ----------------------------------------------------
## Slope estimates of temperature's impact on shape

approach_labels <- c(
  "Ratio"    = "wing / mass",
  "Ratio2"   = "wing² / mass",
  "Resid_ols"= "OLS residuals",
  "Ryding"   = "mass as covariate",
  "SEM_wing" = "SEM (wing)",
  "Sli_est"  = "SLI estimated",
  "Sli_iso"  = "SLI isometry"
)

Num_spp <- length(unique(parms_df$species_))

# Exclude SEM_tarsus from this wing-focused plot
parms_df_p <- parms_df %>%
  filter(term == "B.Tavg", Approach != "SEM_tarsus") %>%
  left_join(Spp_metadata2[, c("species_", "Direction", "Sig_trait")]) %>%
  mutate(Approach = factor(Approach),
         species  = str_replace(species_, "_", "\n")) %>% 
  # Remove large standard errors that make plotting challenging
  filter(std.error < .5)

shape_scale <- c("Both" = 15, "Mass" = 16, "Neither" = 17, "Wing" = 18)

Direction_effect <- unique(Spp_metadata2$Direction)
Direction_effect <- setNames(Direction_effect, Direction_effect)

plot_shape <- function(df, title = NULL, legend = TRUE, drop_y = FALSE) {
  p <- df %>%
    mutate(Sig_trait = str_to_sentence(Sig_trait)) %>%
    ggplot(aes(x = species, y = estimate, color = Approach,
               group = interaction(species, Approach))) +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                  alpha = .8, width = 0,
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(shape = Sig_trait),
               size = 2, position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = NULL, y = expression(beta[T] ~ "on wingyness"),
         title = title) +
    theme(
      axis.text.x = element_text(vjust = .58, angle = 60),
      legend.position = "top"
    ) +
    scale_shape_manual(values = shape_scale) +
    scale_color_discrete(labels = approach_labels)

  if (drop_y) p <- p + theme(axis.title.y = element_blank())
  if (!legend) p <- p + theme(legend.position = "none")
  return(p)
}

Shape_plots <- imap(Direction_effect, \(direction, name) {
  parms_filt <- parms_df_p %>% filter(Direction == direction)
  drop_y     <- direction %in% c("Inverse", "Mixed")
  plot_shape(df = parms_filt, title = name, drop_y = drop_y)
})

# Create common legend
com.leg <- get_legend(plot_shape(df = parms_df_p, title = "All"))

# Plot
ggarrange(
  Shape_plots[[1]], Shape_plots[[2]], Shape_plots[[3]], 
  common.legend = TRUE, legend.grob = com.leg, labels = "auto"
)

# Save
ggsave("Figures/Atlantic_birds_shape.png", bg = "white", height = 7, width = 9)

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

stop()

# Interpretation -------------------------------------------------------------

# This doesn't seem to be a great dataset to test this particular question because there are very few species that follow Bergmann's rule and have a significant correlation between wing and mass. Joliceur 1990 argues that the correlation should be greater than 0.6, which is extremely high. Furthermore, the relationships between temperature and size are pretty small, so temperature may have little effect on size and shape.
#NOTE: I tried with tarsus as well and got similar results

# Next steps --------------------------------------------------------------

# Could try this looking at how shape changes with the amount of forest or the habitat where the bird was captured
# Try with year instead of average temperature. Climate change / deforestation (also increasing through time) may be stronger influences than temperature changes by geography
# Could try with some sort of bill measurement

# EXTRAS ------------------------------------------------------------------
# Test effect of approach -------------------------------------------------

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
# OLS assumptions
ols_mod_l <- map(nj_df_l, \(df){
  ols_mod <- lm(log_wing ~ log_mass + B.Tavg, data = df)
})

# Some departure from homoskedasticity
map(ols_mod_l, \(ols_mod){
  plot(ols_mod, which = 1)
  plot(ols_mod, which = 2)
})

# SMA assumptions
sma_mod_l <- map(nj_df_l, \(df){
  sma_mod <- sma(log_wing ~ log_mass, data = df, method = "SMA")
})

map(sma_mod_l, \(sma_mod){
  plot(sma_mod, which = "residual")
  plot(sma_mod, which = "qq")
})
