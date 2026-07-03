## Relative appendage length key functions script

# Load required libraries
library(MASS)
library(MBESS) # cor2cov function
library(tidyverse)

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
gen_data <- function(n = 3000,
                     b_avg_12 = 0.33,
                     r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                     mean_mass = 80, mean_append = 180, mean_temp = 1,
                     sd_log_morph = 0.10,
                     sd_temp = 0.18,
                     meas_error = 0,
                     transient_error_mass = 0,
                     transient_error_append = 0) {
  
  # Generate log-scale covariance matrix
  Sigma <- gen_cov_mat(b_avg_12 = b_avg_12,
                       r_12 = r_12, r_13 = r_13, r_23 = r_23,
                       sd_log_morph = sd_log_morph,
                       sd_temp = sd_temp)
  
  mu <- c(log(mean_append), log(mean_mass), mean_temp)
  
  sim_log <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma, empirical = TRUE)
  colnames(sim_log) <- c("log_Append", "log_Mass", "Temp")
  
  sim <- as_tibble(sim_log) %>%
    mutate(Append = exp(log_Append),
           Mass = exp(log_Mass),
           Temp_inc = Temp)
  
  # Add error to morphometrics on raw scale
  sd_append <- sd(sim$Append)
  sd_mass <- sd(sim$Mass)
  
  sim2 <- sim %>% rm_outliers(metric = Append, sd_append) %>% 
    rm_outliers(metric = Mass, sd_mass)
  
  sim3 <- sim2 %>%
    mutate(
      Append = Append + rnorm(nrow(sim2), 0, sd = sd_append * meas_error) +
        rnorm(nrow(sim2), 0, sd = sd_append * transient_error_append),
      Mass = Mass + rnorm(nrow(sim2), 0, sd = sd_mass * meas_error) +
        rnorm(nrow(sim2), 0, sd = sd_mass * transient_error_mass),
      Append_log = log(Append),
      Mass_log = log(Mass)
    )
  
  format_temp(sim3)
}

## Generate the variance-covariance matrix (on the log scale) that goes into gen_data function
# vary_sd: vary_sd argument = TRUE allows specification of standard deviation of either mass or appendage from a uniform distribution 
# vary: Given the standard deviations covary together, the vary argument allows the generated sd_log_morph to be assigned to either mass or appendage, and then the unknown standard deviation is solved for algebraically 
gen_cov_mat <- function(b_avg_12 = 0.33,
                        r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                        vary = c("mass", "append"),
                        sd_log_morph = .07,
                        vary_sd = TRUE,
                        sd_temp = 0.18) {
    
    if(vary_sd == TRUE){
      sd_log_morph <- runif(1, 0.05, 0.09)
    }
  
  vary <- match.arg(vary)
  
  # Compute b_OLS
  b_ols <- (2 * b_avg_12 * r_12) / (r_12 + 1)
  
  if (vary == "mass") {
    sd_log_mass <- sd_log_morph
    sd_log_append <- abs(b_ols / r_12 * sd_log_mass)
  } else {
    sd_log_append <- sd_log_morph
    sd_log_mass <- abs(r_12 / b_ols * sd_log_append)
  }
  
  sds <- c(sd_log_append, sd_log_mass, sd_temp)
  
  R <- matrix(c(
    1,     r_12,  r_13,
    r_12,  1,     r_23,
    r_13,  r_23,  1
  ), nrow = 3, byrow = TRUE)
  
  MBESS::cor2cov(cor.mat = R, sd = sds)
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
# Returns one row per group combination (age × sex) with n, b_ols, and p_mw.
# Printed for user inspection: groups with b_ols <= threshold or p_mw > threshold
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
gen_cor_vars <- function(r_12, mu_append, mu_mass, sd_append, sd_mass, transient_error_append, transient_error_mass, meas_error){
  cov_12 <- r_12 * sd_append * sd_mass
  var_cov <- matrix(c(
    sd_append^2, cov_12,
    cov_12, sd_mass^2
  ), nrow = 2)
  
  df <- MASS::mvrnorm(n = 3000, mu = c(mu_append, mu_mass), Sigma = var_cov, empirical = TRUE)
  colnames(df) <- c("Appendage", "Mass")
  
  meas_error_mass <- rnorm(3000, 0, sd = sd_mass * meas_error)
  meas_error_append <- rnorm(3000, 0, sd = sd_append * meas_error)
  # Transient fluctuation (biological error) to mass
  transient_error_mass <- rnorm(3000, 0, sd = sd_mass * transient_error_mass)
  transient_error_append <- rnorm(3000, 0, sd = sd_append * transient_error_append)
  
  df <- as_tibble(df) %>% mutate(
     Mass = Mass + meas_error_mass + transient_error_mass,
     Appendage = Appendage + meas_error_append + transient_error_append
  )
  return(df)
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
