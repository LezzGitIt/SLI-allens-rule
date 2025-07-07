## Key functions script

# Libraries
library(MASS)
library(MBESS) # cor2cov function
library(tidyverse)

# Steps taken to fix: 
# Removed error: .02 and .04
# transient_error in gen_ex_data function
# Switched wing and mass in var-cov matrix: sd_wing <- sd_mass / abs(b_sma_12)
# Removed slice_sample
# Played around with correlations: r_12 <- c(.45, .6, .75), r_13 <- c(-.01, -.025, -.04)
# Simulate on absolute scale, then log transform

# OLD - Function to generate datasets
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
gen_data <- function(
    n = 3000,
    b_avg_12 = 0.33, r_12 = 0.3, 
    means = c(mean_append = 180, mean_mass = 80, mean_temp = 1),
    r_13 = -0.1, r_23 = -0.1,
    meas_error = 0, 
    transient_error_mass = 0, transient_error_append = 0
    ) {
  
  # Covariance matrix 
  var_cov_mat <- gen_3var_cov_abs(b_avg_12, mean_append = means[1], mean_mass = means[2], r_12 = r_12, r_13 = r_13, r_23 = r_23)
  
  raw_data <- MASS::mvrnorm(n = n, mu = means, Sigma = var_cov_mat, empirical = TRUE)
  colnames(raw_data) <- c("Append", "Mass", "Temp_inc")
  raw_data <- as_tibble(raw_data)
  
  # Calc std deviations of mass & append
  sd_mass <- sd(raw_data$Mass)
  sd_append <- sd(raw_data$Append)
  
  ## Add various sources of error
  # Add measurement error to mass and append (2% of std deviation)
  meas_error_mass <- rnorm(3000, 0, sd = sd_mass * meas_error)
  meas_error_append <- rnorm(3000, 0, sd = sd_append * meas_error)
  # Transient fluctuation (biological error) to mass
  transient_error_mass <- rnorm(3000, 0, sd = sd_mass * transient_error_mass)
  transient_error_append <- rnorm(3000, 0, sd = sd_append * transient_error_append)
  
  raw_data <- raw_data %>%
    mutate(Mass = Mass + meas_error_mass + transient_error_mass,
           Append = Append + meas_error_append + transient_error_append,
           Append_log = log(Append),
           Mass_log = log(Mass),
           Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)) %>%
    arrange(Temp_inc) %>%
    #slice_sample(n = 200, by = Temp_bin) %>% 
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:4) ~ 4,  
      Temp_bin %in% c(12:15) ~ 12,
      .default = Temp_bin
    ))
  if(any(raw_data$Mass < 0) | any(raw_data$Append < 0)){
    print("Removing morphometrics < 0")
    raw_data %>% filter(Mass > 0 & Append > 0)
  }
  return(raw_data)
}


# Simulate data for a species with isometric scaling (log-log slope = 1)
library(tidyverse)
library(MASS)
library(lmodel2)

# Simulate data
dat <- gen_data(b_avg_12 = 1) # Isometry

# Check summary stats: raw and log SDs
dat %>% summarize(
  mean_mass = mean(Mass),
  mean_append = mean(Append),
  sd_mass_raw = sd(Mass),
  sd_append_raw = sd(Append),
  sd_mass_log = sd(log(Mass)),
  sd_append_log = sd(log(Append)),
  sd_ratio_log = sd(log(Append)) / sd(log(Mass))
)

# Estimate OLS slope on log-log scale
ols <- lm(Append_log ~ Mass_log, data = dat)
summary(ols)$coef

# Estimate SMA slope
sma <- lmodel2(log(Append) ~ log(Mass), data = dat)
sma$regression.results



if(FALSE){
  nightjar <- gen_data(b_avg_12 = 0.303)
  whippoor <- gen_data(b_sma_12 = 0.4)
  nighthawk <- gen_data(b_sma_12 = 0.377)
  
  map(list(nightjar, whippoor, nighthawk), 
      ~ with(.x, sd(Wing_log) / sd(Mass_log))) 
}

format_temp <- function(df){
  df %>% mutate(
    Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)) %>%
    arrange(Temp_inc) %>%
    #slice_sample(n = 200, by = Temp_bin) %>% 
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:4) ~ 4,  
      Temp_bin %in% c(12:15) ~ 12,
      .default = Temp_bin
    ))
}



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
  
  # Add error on raw scale
  sd_append <- sd(sim$Append)
  sd_mass <- sd(sim$Mass)
  
  # Morphometrics 
  sim <- sim %>%
    mutate(
      Append = Append + rnorm(n, 0, sd = sd_append * meas_error) +
        rnorm(n, 0, sd = sd_append * transient_error_append),
      Mass = Mass + rnorm(n, 0, sd = sd_mass * meas_error) +
        rnorm(n, 0, sd = sd_mass * transient_error_mass),
      Append_log = log(Append),
      Mass_log = log(Mass)
    )
  
  format_temp(sim)
}



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



# OLD - Generate data on the raw (absolute) scale 
gen_3var_cov_abs <- function(b_sma_12 = 0.33, 
                             mean_mass = 80, 
                             mean_wing = 180, 
                             mean_temp = 1, 
                             sd_log_mass = 0.10, 
                             sd_temp = 0.2,
                             r_12 = 0.3, 
                             r_13 = -0.1, 
                             r_23 = -0.1) {
  
  # The relationship between b_sma & the SDs of wing and mass is deterministic. Calculate sd_wing (log scale) from b_sma and sd_mass
  sd_log_wing <- abs(b_sma_12) * sd_log_mass
  
  # Convert log-scale SDs to raw-scale SDs for lognormal
  var_mass <- (exp(sd_log_mass^2) - 1) * exp(2 * log(mean_mass) + sd_log_mass^2)
  var_wing <- (exp(sd_log_wing^2) - 1) * exp(2 * log(mean_wing) + sd_log_wing^2)
  
  sd_mass <- sqrt(var_mass)
  sd_wing <- sqrt(var_wing)
  
  # Examine cor2cov function
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

# In progress
gen_3var_cov_abs <- function(b_avg_12 = 0.33, 
                             mean_append = 180,
                             mean_mass = 80, 
                             sd_log_mass = 0.10, 
                             sd_temp = 0.2,
                             r_12 = 0.3, 
                             r_13 = -0.3, 
                             r_23 = -0.3) {
  
  # Ratio of standard deviations that gives desired average of B_ols and B_sma
  ratio_append_mass <- 2 * b_avg_12 / (r_12 + 1 / r_12)
  
  # Set log-scale SDs for lognormal
  sd_log_append <- ratio_append_mass * sd_log_mass
  
  # Convert log-scale SDs to raw-scale SDs for lognormal
  var_mass <- (exp(sd_log_mass^2) - 1) * exp(2 * log(mean_mass) + sd_log_mass^2)
  var_append <- (exp(sd_log_append^2) - 1) * exp(2 * log(mean_append) + sd_log_append^2)
  
  sd_mass <- sqrt(var_mass)
  sd_append <- sqrt(var_append)
  
  # Covariances
  cov_12 <- r_12 * sd_append * sd_mass
  cov_13 <- r_13 * sd_append * sd_temp
  cov_23 <- r_23 * sd_mass * sd_temp
  
  Sigma <- matrix(c(
    sd_append^2, cov_12, cov_13,
    cov_12,      sd_mass^2, cov_23,
    cov_13,      cov_23,    sd_temp^2
  ), nrow = 3)
  
  if (any(eigen(Sigma, only.values = TRUE)$values <= 0)) {
    warning("Covariance matrix is not positive definite")
  }
  return(Sigma)
}

gen_cov_mat <- function(b_avg_12 = 0.33,
                        r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                        sd_log_mass = 0.07,
                        vary_sd = TRUE,
                        sd_temp = 0.20) {
  
  if(vary_sd == TRUE){
    sd_log_mass <- runif(1, 0.05, 0.09)
  }
  
  # Solve for log-scale SD of appendage
  b_ols <- (2 * b_avg_12 * r_12) / (r_12 + 1)
  sd_log_append <- abs(b_ols / r_12 * sd_log_mass)
  
  # Define correlation matrix
  R <- matrix(c(
    1,     r_12,  r_13,
    r_12,  1,     r_23,
    r_13,  r_23,  1
  ), nrow = 3, byrow = TRUE)
  
  # Corresponding SDs for log(Append), log(Mass), log(Temp)
  sds <- c(sd_log_append, sd_log_mass, sd_temp)
  
  # Compute covariance matrix
  Sigma <- MBESS::cor2cov(cor.mat = R, sd = sds)
  print(sds[1])
  return(Sigma)
}

# Testing
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


# SLI function
calc_sli <- function(df, b_sli = 0.33, rename_col = FALSE){
  # L0 is the average mass, essentially allowing for comparison of wing lengths for a given mass
  L0 <- mean(df$Mass)
  df_sli <- df %>% mutate(sli = Append * (L0 / Mass)^b_sli) %>%
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

# Switch to cor2cov?
cor_mat <- matrix(c(1, r_12, r_12, 1), nrow = 2)

#%>% #rename_with(.cols = .,.fn = ~c("Wing", "Mass"))
