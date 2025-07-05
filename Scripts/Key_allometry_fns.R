## Key functions script

# Libraries
library(MASS)
library(tidyverse)

# Steps taken to fix: 
# Removed error: .02 and .04
# transient_error in gen_ex_data function
# Switched wing and mass in var-cov matrix: sd_wing <- sd_mass / abs(b_sma_12)
# Removed slice_sample
# Played around with correlations: r_12 <- c(.45, .6, .75), r_13 <- c(-.01, -.025, -.04)
# Simulate on absolute scale, then log transform

# Function to generate datasets
gen_data <- function(b_sma_12 = NULL, r_12, r_13, r_23, meas_error = 0.00, transient_error = 0.00) {
  Sigma <- gen_3var_cov(b_sma_12, r_12, r_13, r_23)
  data <- mvrnorm(n = 3000, mu = c(45, 15, 1.3), Sigma = Sigma, empirical = TRUE)
  morph_temp <- tibble(Wing = data[,1], # Wing = appendage
                       Mass= data[,2], # Mass = body size
                       Temp_inc= data[,3])
  
  # Calc std deviations of mass & wing
  sd_mass <- sd(morph_temp$Mass)
  sd_wing <- sd(morph_temp$Wing)
  
  ## Add various sources of error
  # Add measurement error to mass and wing (2% of std deviation)
  meas_error_mass <- rnorm(3000, 0, sd = sd_mass * meas_error)
  meas_error_wing <- rnorm(3000, 0, sd = sd_wing * meas_error)
  # Transient fluctuation (biological error) to mass
  transient_error_mass <- rnorm(3000, 0, sd = sd_mass * transient_error)
  
  morph_temp2 <- morph_temp %>% 
    mutate(Mass = Mass + meas_error_mass + transient_error_mass,
           Wing = Wing + meas_error_wing) %>% 
    mutate(Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, 
                          ordered_result = TRUE)) %>%
    arrange(Temp_inc) %>%
    slice_sample(n = 200, by = Temp_bin) %>% 
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:4) ~ 4,  
      Temp_bin %in% c(12:15) ~ 12,
      .default = Temp_bin
    ))
  return(morph_temp2)
}



# In progress gen_data 
gen_data <- function(b_sma_12 = 0.33, 
                     r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                     n = 3000) {
  
  cov_result <- gen_3var_cov(b_sma_12, r_12 = r_12, r_13 = r_13, r_23 = r_23)
  Sigma <- cov_result$Sigma
  means <- cov_result$means
  
  raw_data <- MASS::mvrnorm(n = n, mu = means, Sigma = Sigma)
  colnames(raw_data) <- c("Wing", "Mass", "Temp_inc")
  
  as_tibble(raw_data) %>%
    mutate(Wing_log = log(Wing),
           Mass_log = log(Mass),
           Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)) %>%
    arrange(Temp_inc) %>%
    slice_sample(n = 200, by = Temp_bin) %>% 
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:4) ~ 4,  
      Temp_bin %in% c(12:15) ~ 12,
      .default = Temp_bin
    ))
}

nightjar <- gen_data(b_sma_12 = 0.303)
whippoor <- gen_data(b_sma_12 = 0.4)
nighthawk <- gen_data(b_sma_12 = 0.377)

map(list(nightjar, whippoor, nighthawk), 
    ~ with(.x, sd(log_Wing) / sd(log_Mass)))


# OLD copy of working , but with variance far too high 
# Generate covariance matrix 
gen_3var_cov <- function(b_sma_12 = NULL, r_12 = .5, r_13 = 0, r_23 = 0) { 
  if(!is.null(b_sma_12)) { 
    S1 <- abs(b_sma_12) 
  } else { 
    S1 <- 1
  } 
  S2 <- 1
  cov_12 <- r_12 * S1 * S2 
  cov_13 <- r_13 * S1 * 1 
  cov_23 <- r_23 * S2 * 1 
  matrix(c(S1^2, cov_12, cov_13, cov_12, S2^2, cov_23, cov_13, cov_23, 1^2), nrow = 3) 
}


#if(FALSE){
# In progress 1 -- generate covariance matrix 
gen_3var_cov <- function(b_sma_12 = NULL, r_12 = .5, r_13 = 0, r_23 = 0) { 
  var_factor <- 0.8
  if(!is.null(b_sma_12)) { 
    S1 <- abs(b_sma_12) 
  } else { 
    S1 <- var_factor
  } 
  S2 <- var_factor
  cov_12 <- r_12 * S1 * S2 #  S2 was 1 previously 
  cov_13 <- r_13 * S1 * 1 
  cov_23 <- r_23 * S2 * 1 
  matrix(c(S1^2, cov_12, cov_13, cov_12, S2^2, cov_23, cov_13, cov_23, var_factor^2), nrow = 3) 
}

# In progress 2 -- generate covariance matrix 
gen_3var_cov <- function(b_sma_12 = NULL, 
                         sd_wing = 4, 
                         sd_temp = 0.3,
                         r_12 = 0.5, 
                         r_13 = 0, 
                         r_23 = 0) {
  
  sd_mass <- sd_wing / abs(b_sma_12)
  
  b_sma = sd_wing / sd_mass
  1/b_sma = sd_mass / sd_wing
  
  cov_12 <- r_12 * sd_wing * sd_mass
  cov_13 <- r_13 * sd_wing * sd_temp
  cov_23 <- r_23 * sd_mass * sd_temp
  
  Sigma <- matrix(c(
    sd_wing^2, cov_12, cov_13,
    cov_12,   sd_mass^2, cov_23,
    cov_13,   cov_23,    sd_temp^2
  ), nrow = 3)
  
  if (any(eigen(Sigma, only.values = TRUE)$values <= 0)) {
    warning("Covariance matrix is not positive definite")
  }
  
  return(Sigma)
}
#}

# Raw scale 
gen_3var_cov_abs <- function(b_sma_12 = 0.33, 
                             mean_mass = 80, 
                             mean_wing = 180, 
                             mean_temp = 1.3, 
                             sd_log_mass = 0.10, 
                             sd_temp = 0.2,
                             r_12 = 0.3, 
                             r_13 = -0.1, 
                             r_23 = -0.1) {
  
  # Set sd on log scale
  sd_log_wing <- b_sma_12 * sd_log_mass
  
  # Convert log-scale SDs to raw-scale SDs for lognormal
  var_mass <- (exp(sd_log_mass^2) - 1) * exp(2 * log(mean_mass) + sd_log_mass^2)
  var_wing <- (exp(sd_log_wing^2) - 1) * exp(2 * log(mean_wing) + sd_log_wing^2)
  
  sd_mass <- sqrt(var_mass)
  sd_wing <- sqrt(var_wing)
  
  cov_12 <- r_12 * sd_mass * sd_wing
  cov_13 <- r_13 * sd_wing * sd_temp
  cov_23 <- r_23 * sd_mass * sd_temp
  
  Sigma <- matrix(c(
    sd_wing^2, cov_12, cov_13,
    cov_12,   sd_mass^2, cov_23,
    cov_13,   cov_23,    sd_temp^2
  ), nrow = 3)
  
  if (any(eigen(Sigma, only.values = TRUE)$values <= 0)) {
    warning("Covariance matrix is not positive definite")
  }
  
  list(Sigma = Sigma, means = c(mean_wing, mean_mass, mean_temp))
}


# SLI function
calc_sli <- function(df, b_sli = 0.33, rename_col = FALSE){
  # L0 is the average mass, essentially allowing for comparison of wing lengths for a given mass
  L0 <- mean(df$Mass)
  df_sli <- df %>% mutate(sli = Wing * (L0 / Mass)^b_sli) %>%
    arrange(desc(sli))
  if(rename_col != FALSE){df_sli <- df_sli %>% rename( {{ rename_col }} := sli)}
  return(df_sli)
}

# calc_lambda function
calc_lambda <- function(x, y){ 
  cv_y <- sd({{ y }}) / mean({{ y }}) 
  cv_x <- sd({{ x }}) / mean({{ x }}) 
  (cv_y^2) / (cv_x^2) 
  } 