## Relative appendage length key functions script

### Migration to the sliR package is in progress. Functions ported to sliR keep their original names, arguments and output columns here, so no call site changes; their bodies are thin wrappers. Everything not in sliR (format_temp, rm_outliers, run_sma_mod, format_sma_parms, gen_ex_data, calc_lambda, classify_direction) stays a local helper.
# Migrated so far: gen_data() [-> sliR::sim_allometric()], gen_cov_mat() [-> sliR::build_cov_mat()], gen_cor_vars() [-> sliR::sim_correlated()]
# Still to migrate: calc_sli(), build_sli_slopes_tbl(), build_group_cor_tbl()
# remotes::install_github("LezzGitIt/sliR")

# Load required libraries
# MASS is no longer used by this file, but is left attached because supplementary_info.qmd sources this script without loading MASS itself; dropping it here would change that document's search path.
library(MASS)
library(tidyverse)
library(sliR)

# Creation of temperature bins for plotting and examination of scaling intercepts
format_temp <- function(df){
  df %>% mutate(
    Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)) %>%
    arrange(Temp_inc) %>%
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:5) ~ 5,  
      Temp_bin %in% c(11:15) ~ 11,
      .default = Temp_bin
    )) #%>% slice_sample(n = 200, by = Temp_bin) 
}
?cut_number # Consider cut_number to make groups with equal number of individuals

# Remove outliers 3 or more SDs from mean
rm_outliers <- function(df, metric, sd_metric) {
  df %>% filter(!({{ metric }} > mean({{ metric }}) + 3 * {{ sd_metric }} |
                    {{ metric }} < mean({{ metric }}) - 3 * {{ sd_metric }}))
}

## Generate data (on the log-scale) using var-cov matrix
# Specify measurement error and transient 'error'
### Now a thin wrapper over sliR::sim_allometric(). Arguments and output columns are unchanged, so every call site still works. Temperature is drawn as a normal gradient, matching the multivariate-normal draw this function used to do; sliR's own default gradient is uniform.
# vary_sd was previously hidden inside gen_cov_mat(), where it silently overwrote whatever sd_log_morph the caller passed. It is surfaced here as an argument, defaulting to TRUE to preserve the historical behaviour. Setting vary_sd = FALSE honours sd_log_morph. It affects raw-scale dispersion and skew only: the log-scale slopes and correlations are invariant to it.
gen_data <- function(n = 3000,
                     b_avg_12 = 0.33,
                     r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                     mean_mass = 80, mean_append = 180, mean_temp = 1,
                     sd_log_morph = 0.10,
                     sd_temp = 0.18,
                     vary_sd = TRUE,
                     meas_error = 0,
                     transient_error_mass = 0,
                     transient_error_append = 0) {

  if (vary_sd) sd_log_morph <- runif(1, 0.05, 0.09)

  sim <- sliR::sim_allometric(
    n             = n,
    b_avg         = b_avg_12,
    r_app_mass    = r_12,
    sd_log_mass   = sd_log_morph,
    mean_append   = mean_append,
    mean_mass     = mean_mass,
    gradient      = "Temp_inc",
    gradient_dist = "normal",
    mean_gradient = mean_temp,
    sd_gradient   = sd_temp,
    r_grad_app    = r_13,
    r_grad_mass   = r_23,
    meas_error             = meas_error,
    transient_error_append = transient_error_append,
    transient_error_mass   = transient_error_mass,
    trim_sd = 3
  )

  format_temp(sim)
}

## Generate the variance-covariance matrix (on the log scale), for log(appendage), log(mass) and temperature
### Now a thin wrapper over sliR::build_cov_mat(). The allometry has two degrees of freedom, so b_avg_12 and r_12 together pin both slopes, and one standard deviation then fixes the scale; sliR solves for the other. Output is an unnamed 3x3 matrix, as before.
# vary: which morphological trait sd_log_morph refers to. The two standard deviations covary, so naming one solves for the other.
# vary_sd: TRUE draws sd_log_morph from Uniform(0.05, 0.09) and ignores whatever was passed. Kept for backwards compatibility. It affects raw-scale dispersion and skew only: the log-scale slopes and correlations are invariant to it.
gen_cov_mat <- function(b_avg_12 = 0.33,
                        r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                        vary = c("mass", "append"),
                        sd_log_morph = .07,
                        vary_sd = TRUE,
                        sd_temp = 0.18) {

  vary <- match.arg(vary)
  if (vary_sd) sd_log_morph <- runif(1, 0.05, 0.09)

  sd_arg <- if (vary == "mass") list(sd_log_mass = sd_log_morph) else list(sd_log_append = sd_log_morph)

  Sigma <- do.call(sliR::build_cov_mat, c(
    list(b_avg = b_avg_12, r_app_mass = r_12,
         gradient = "Temp", r_grad_app = r_13, r_grad_mass = r_23),
    sd_arg
  ))

  ## build_cov_mat() always places the gradient at unit SD, because it describes a correlation structure rather than temperature's units. Rescale that block to sd_temp.
  scale_temp <- diag(c(1, 1, sd_temp))
  unname(scale_temp %*% Sigma %*% scale_temp)
}


# Fit an SMA model of Append_log ~ Mass_log, optionally allowing the slope to
# vary with binned temperature (Temp_bin). Used to validate simulated species'
# direction of shape change via their SMA intercepts.
run_sma_mod <- function(df, interaction = FALSE) {
  if (!interaction) {
    smatr::sma(Append_log ~ Mass_log + Temp_bin, data = df, method = "SMA")
  } else {
    smatr::sma(Append_log ~ Mass_log * Temp_bin, data = df, method = "SMA")
  }
}

# Tidy the per-Temp_bin coefficients (intercept/slope) from an sma() model
# fit with run_sma_mod(), decoding the bin label back to a numeric/label pair.
format_sma_parms <- function(sma_mod) {
  coef(sma_mod) %>%
    tibble::rownames_to_column("Temp_inc") %>%
    dplyr::mutate(
      Temp_inc = stringr::str_pad(Temp_inc, side = "left", width = 2, pad = "0"),
      Temp_inc = stringr::str_replace(Temp_inc, "^([0-9])([0-9])$", "\\1.\\2"),
      Temp_label = paste0(Temp_inc, "°C"),
      Temp_inc = as.numeric(Temp_inc)
    ) %>%
    tibble::tibble()
}

# Regenerate raw individual-level data for one or more hypothetical species
# (rows of a Parms_mat-style tibble) via gen_data(), for illustrative figures.
gen_ex_data <- function(Parms_mat, transient_error_mass = 0, transient_error_append = 0) {
  Cols <- Parms_mat %>% dplyr::select(dplyr::starts_with(c("b_", "r_")))
  Parms_mat %>%
    dplyr::mutate(coefs = purrr::pmap(Cols, \(...) gen_data(...,
                                              transient_error_mass   = transient_error_mass,
                                              transient_error_append = transient_error_append))) %>%
    tidyr::unnest(coefs)
}

# Build a summary tibble of per-group SMA slopes for SLI estimation.
# Fits one SMA model per element of `control` (e.g. Age, Sex separately),
# then averages the resulting slopes for every group combination.
# Rows with unknown-coded values or NA in control variables are excluded from
# slope fitting but all combinations present in the known data are represented.
# Mass argument is NSE (default = Mass) for compatibility with lowercase column names.
# Requires smatr and rlang to be loaded (done by the calling script).
build_sli_slopes_tbl <- function(df, Append, Mass = Mass, control,
                                  unknown_codes = c("Unk", "U", "Unknown")) {
  app_nm  <- rlang::as_label(rlang::enquo(Append))
  mass_nm <- rlang::as_label(rlang::enquo(Mass))

  df_fit <- df %>%
    mutate(.log_app = log(.data[[app_nm]]), .log_mass = log(.data[[mass_nm]])) %>%
    filter(if_all(all_of(control), \(x) !is.na(x) & !x %in% unknown_codes))

  # One SMA per control variable; collect per-level slopes
  slope_tbls <- map(control, \(grp_var) {
    fmla <- as.formula(paste(".log_app ~ .log_mass *", grp_var))
    fit  <- smatr::sma(fmla, data = df_fit, method = "SMA")
    as_tibble(coef(fit), rownames = grp_var) %>%
      rename(!!paste0("b_sma_", grp_var) := slope) %>%
      dplyr::select(all_of(grp_var), starts_with("b_sma_"))
  })

  # All group combinations present in the known data, with sample sizes
  n_tbl <- df_fit %>% count(!!!syms(control))

  # Join per-variable slopes onto the combination table; average for b_sli_avg
  slopes_tbl <- n_tbl
  for (i in seq_along(control)) {
    slopes_tbl <- left_join(slopes_tbl, slope_tbls[[i]], by = control[[i]])
  }
  slope_cols <- paste0("b_sma_", control)
  slopes_tbl %>%
    mutate(b_sli_avg = rowMeans(across(all_of(slope_cols)), na.rm = FALSE))
}

# Per-group mass ~ appendage OLS correlation table.
# Returns one row per group combination (age × sex) with n, b_ols, r_mw, and p_mw.
# Printed for user inspection: groups with r_mw <= cor_min or p_mw > threshold
# lack a meaningful allometric relationship and should not drive per-group SMA slopes.
build_group_cor_tbl <- function(df, Append, Mass = Mass, control,
                                 unknown_codes = c("Unk", "U", "Unknown")) {
  app_nm  <- rlang::as_label(rlang::enquo(Append))
  mass_nm <- rlang::as_label(rlang::enquo(Mass))

  df_fit <- df %>%
    filter(if_all(all_of(control), \(x) !is.na(x) & !x %in% unknown_codes))

  df_fit %>%
    group_by(!!!syms(control)) %>%
    group_modify(\(grp, key) {
      fmla <- as.formula(paste(mass_nm, "~", app_nm))
      m    <- lm(fmla, data = grp)
      tibble(
        n     = nrow(grp),
        b_ols = coef(m)[[app_nm]],
        r_mw  = cor(grp[[app_nm]], grp[[mass_nm]], use = "complete.obs"),
        p_mw  = summary(m)$coefficients[app_nm, "Pr(>|t|)"]
      )
    })
}

# calc_sli function: see Peig & Green (2009)
# When control = NULL: scalar b_sli applied to all rows.
# When control is a character vector (e.g. c("Age","Sex")): build_sli_slopes_tbl()
# estimates per-group slopes and each individual gets their group's averaged slope.
# Unknown-coded rows (in control variables) receive sli = NA.
# Mass argument is NSE (default = Mass) for compatibility with lowercase column names.
calc_sli <- function(df, Append = Append, Mass = Mass, b_sli = 0.33,
                     rename_col = FALSE, control = NULL) {
  mass_nm <- rlang::as_label(rlang::enquo(Mass))
  L0      <- mean(df[[mass_nm]], na.rm = TRUE)

  if (!is.null(control)) {
    app_q  <- rlang::enquo(Append)
    mass_q <- rlang::enquo(Mass)
    slopes_tbl <- build_sli_slopes_tbl(df, Append = !!app_q, Mass = !!mass_q,
                                        control = control)
    df_sli <- df %>%
      left_join(slopes_tbl %>% dplyr::select(all_of(control), b_sli_avg), by = control) %>%
      mutate(sli = {{ Append }} * (L0 / {{ Mass }})^b_sli_avg) %>%
      dplyr::select(-b_sli_avg)
  } else {
    df_sli <- df %>%
      mutate(sli = {{ Append }} * (L0 / {{ Mass }})^b_sli)
  }

  df_sli <- df_sli %>% arrange(desc(sli))
  if (rename_col != FALSE) df_sli <- df_sli %>% rename({{ rename_col }} := sli)
  return(df_sli)
}

# calc_lambda function: calculate the empirical coefficients of variation 
calc_lambda <- function(x, y){ 
  cv_y <- sd({{ y }}) / mean({{ y }}) 
  cv_x <- sd({{ x }}) / mean({{ x }}) 
  (cv_y^2) / (cv_x^2) 
} 

## Generate correlated mass and wing data
# Helpful for playing around to see how slopes vary with different relationships of mass and wing
# var1 = appendage, mass = mass
### Now a thin wrapper over sliR::sim_correlated(). Raw scale, no log transform and no allometric target, so it is the counterpart to gen_data() for building intuition about how a fitted slope responds to correlation and to error on one trait but not the other.
# sliR names the appendage column Append; rename it back to Appendage so existing call sites are unaffected. n was hardcoded at 3000 and is now an argument, defaulting to that.
gen_cor_vars <- function(r_12, mu_append, mu_mass, sd_append, sd_mass,
                         transient_error_append, transient_error_mass, meas_error,
                         n = 3000) {
  sliR::sim_correlated(
    n         = n,
    r         = r_12,
    mu_append = mu_append,
    mu_mass   = mu_mass,
    sd_append = sd_append,
    sd_mass   = sd_mass,
    meas_error             = meas_error,
    transient_error_append = transient_error_append,
    transient_error_mass   = transient_error_mass
  ) %>%
    dplyr::rename(Appendage = Append)
}

# Classify shapeshifting direction (Bergmann's, Inverse Bergmann's, Mixed, Stable)
# from tidy lm output with mass and wing gradient models per species.
# mass_dir = TRUE means mass decreases along gradient (Bergmann's direction for mass).
# wing_dir = TRUE means wing increases along gradient (Allen's direction for wing).
classify_direction <- function(mods_tbl, p_threshold = 0.05,
                               species_col = "species_",
                               mass_dv = "mass", wing_dv = "wing") {
  mods_tbl %>%
    rename(species_ = !!sym(species_col)) %>%
    mutate(sig = p.value < p_threshold) %>%
    group_by(species_) %>%
    summarise(
      n_sig     = sum(sig),
      mass_dir  = estimate[dv == mass_dv] < 0,
      wing_dir  = estimate[dv == wing_dv] > 0,
      mass_sig  = sig[dv == mass_dv],
      wing_sig  = sig[dv == wing_dv],
      Sig_trait = case_when(
        n_sig == 0 ~ "Neither",
        n_sig == 2 ~ "both",
        TRUE       ~ dv[sig]
      ),
      .groups = "drop"
    ) %>%
    mutate(
      sole_decr = case_when(
        mass_sig & !wing_sig ~  mass_dir,
        !mass_sig & wing_sig ~ !wing_dir,
        .default = NA
      ),
      Direction = case_when(
        n_sig == 0                             ~ "Stable",
        n_sig == 1 & sole_decr == TRUE         ~ "Bergmann's",
        n_sig == 1 & sole_decr == FALSE        ~ "Inverse Bergmann's",
        n_sig == 2 &  mass_dir & !wing_dir     ~ "Bergmann's",
        n_sig == 2 & !mass_dir &  wing_dir     ~ "Inverse Bergmann's",
        n_sig == 2 &  mass_dir &  wing_dir     ~ "Mixed - Wingier",
        n_sig == 2 & !mass_dir & !wing_dir     ~ "Mixed - Fatter",
        TRUE ~ "Check"
      )
    ) %>%
    dplyr::select(-sole_decr) %>%
    rename(!!sym(species_col) := species_)
}
