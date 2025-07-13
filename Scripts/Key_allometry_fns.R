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


# calc_sli function: see Peig & Green (2009)
calc_sli <- function(df, b_sli = 0.33, rename_col = FALSE){
  # L0 is the average mass, essentially allowing for comparison of wing lengths for a given mass
  L0 <- mean(df$Mass)
  df_sli <- df %>% mutate(sli = Append * (L0 / Mass)^b_sli) %>%
    arrange(desc(sli))
  if(rename_col != FALSE){df_sli <- df_sli %>% rename( {{ rename_col }} := sli)}
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
