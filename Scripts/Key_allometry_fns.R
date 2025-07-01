## Key functions script

# Libraries
library(MASS)
library(tidyverse)

# Function to generate datasets
gen_data <- function(b_sma_12 = NULL, r_12, r_13, r_23, meas_error = 0.02, transient_error = 0.04) {
  Sigma <- gen_3var_cov(b_sma_12, r_12, r_13, r_23)
  data <- mvrnorm(n = 3000, mu = c(3, 3, 3), Sigma = Sigma, empirical = TRUE)
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
  matrix(c(S1^2, cov_12, cov_13, cov_12, S2^2, cov_23, cov_13, cov_23, 1 ), nrow = 3) 
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
  cv_x <- sd({{ x }}) / mean({{ x }}) 
  cv_y <- sd({{ y }}) / mean({{ y }}) 
  (cv_y^2) / (cv_x^2) } 