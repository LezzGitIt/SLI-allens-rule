## Simulation code to better understand how allometry influences metrics of size & shape shifting
# Goal is to show interesting cases where allometry is important to understanding the true morphological response. Try to... 
# Go beyond just wing & mass (& even beyond just birds)

# Libraries ---------------------------------------------------------------
## Load libraries
library(tidyverse)
library(janitor)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(smatr)
library(broom.mixed)
library(MASS) # mvrnorm() function
library(ggpmisc) # Plots MA or SMA lines of best fit
library(scales)
ggplot2::theme_set(theme_cowplot())
#library(conflicted)
#conflicts_prefer(dplyr::select)
#conflicts_prefer(dplyr::filter)

# Background --------------------------------------------------------------

# b_ols = cov[x,y] / var[x] OR equivalently 
# b_ols = r * Sy / Sx where S is the std. dev
# b_sma = sign(r) * Sy / Sx , ie if r is negative then b_sma is also negative

# NOTE:: From the formulas above, we have 
# b_sma = (r * (Sy / Sx)) / r # rs cancel out


# Simulate functions --------------------------------------------------------
## Generate functions to create datasets with different amounts of correlation (r), OLS slopes, and SMA slopes, and plot 3x3 grid.

## gen_parm_combos() functions
# Goal is to provide either b_ols_vals OR r_vals, & generate the missing values for a given set of b_sma_vals. 
# Formulas used are as follow:
# b_sma = b_ols / r
# b_ols = b_sma * r 
# r = b_ols / b_sma
gen_parm_combos <- function(
    b_sma_vals = c(0.2, 0.33, 0.5), 
    b_ols_vals = NULL, r_vals = NULL) {
  stopifnot(xor(is.null(b_ols_vals), is.null(r_vals)))  # only one should be provided
  
  if (!is.null(r_vals)) {
    parm_combos <- expand.grid(b_sma = b_sma_vals, r = r_vals)
    parm_combos$b_ols <- parm_combos$b_sma * parm_combos$r
  } else {
    parm_combos <- expand.grid(b_sma = b_sma_vals, b_ols = b_ols_vals)
    parm_combos$r <-  parm_combos$b_ols / parm_combos$b_sma
  }
  parm_combos %>% mutate(across(everything(.), ~round(.x, 2)))
}
parm_combos <- gen_parm_combos(r_vals = c(.2, .4, .6)) #b_ols_vals = c(.08, .15, .19)

## simulate_SMA_data() function
# Simulate data under the constraints in parm_combos, & estimate the parameters from the simulated data 
simulate_SMA_data <- function(parm_combos, n = 1000, seed = NULL){
  set.seed(seed)
  sim_list <- list()
  for (i in seq_len(nrow(parm_combos))) {
    print(i)
    b_ols <- parm_combos$b_ols[i]
    r <- parm_combos$r[i]
    b_sma <- parm_combos$b_sma[i]
    
    X <- rnorm(n, mean = 0, sd = 1)
    var_X <- var(X)
    
    # Compute the appropriate amount of noise to generate Y, resulting in the appropriate correlation between X & Y
    var_epsilon <- (b_ols^2 * var_X) * ((1 - r^2) / r^2)
    sd_epsilon <- sqrt(var_epsilon)
    
    eps <- rnorm(n, mean = 0, sd = sd_epsilon)
    Y <- b_ols * X + eps
    
    sim_list[[paste0("b_sma=", b_sma, "_r=", round(r, 2))]] <- tibble(
      X = X,
      Y = Y,
      true_r = r,
      true_b_ols = b_ols,
      true_b_sma = b_sma,
      est_r = cor(X, Y),
      est_b_ols = coef(lm(Y ~ X))[2],
      est_b_sma = coef(sma(Y ~ X))[2]
    )
  }
  bind_rows(sim_list, .id = "condition")
}

# Generate data frame
sim_df <- simulate_SMA_data(parm_combos = parm_combos) %>% 
  mutate(Scaling = case_when(
    true_b_sma == .2 ~ "Hypoallometric",
    true_b_sma == .33 ~ "Isometric",
    true_b_sma == .5 ~ "Hyperallometric",
  ))

# Ensure that estimated parameters are similar to the true parameter values
create_SMA_summary <- function(sim_df) {
  sim_df %>% dplyr::select(starts_with(c("true", "est"))) %>%
    distinct()
}
create_SMA_summary(sim_df)

# Plot 3x3 grid 
# Row 1: Hypoallometry, row 2: isometry, row 3: hyperallometry
# NOTE:: At lower slopes the difference between OLS & SMA regression is greater
plot_SMA_grid <- function(sim_data) {
  sim_data %>%
    group_by(condition) %>%
    mutate(xbar = mean(X), ybar = mean(Y)) %>%
    ungroup() %>%
    ggplot(aes(X, Y)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red") +  # OLS line
    geom_abline(data = sim_data %>% distinct(condition, xbar = mean(X), ybar = mean(Y), est_b_sma), aes(slope = est_b_sma, intercept = ybar - est_b_sma * xbar), color = "blue", size = 1) +
    facet_wrap(~condition, scales = "free") +
    theme_minimal() +
    labs(title = "SMA Regression: Varying Slopes and Correlations") +
    ylim(c(-1, 1))
}
plot_SMA_grid(sim_data = sim_df)

# Simulate 3 vars --------------------------------------------------------
## GOAL:: I want to include a third variable, temperature. So we'll have mass (X1), temperature (X2), and Wing (Y). 
## Example analysis: If we are interested in understanding how temp increase influence wingyness.. We could use 1) residuals from allometric scaling with mass on x-axis (wing_resid ~ temp), or 2) multiple regression (Wing ~ temp + mass). 
# I want to be able to set different correlations between X1 and X2, as well as different b_sma slopes for X1 & Y. 

#gen_parm_combos_mv() function to generate combination of true parameter values
gen_parm_combos_mv <- function(
    b_sma_vals = c(0.2, 0.33, 0.5), 
    r_x1y_vals = NULL, 
    b1_ols_vals = NULL, 
    r12_vals = c(-0.15, -0.3)
) {
  if (is.null(r_x1y_vals) && is.null(b1_ols_vals)) {
    stop("Must supply either r_x1y_vals or b1_ols_vals.")
  }
  if (!is.null(r_x1y_vals) && !is.null(b1_ols_vals)) {
    stop("Supply only one of r_x1y_vals or b1_ols_vals.")
  }
  
  if (!is.null(r_x1y_vals)) {
    parm_combos <- expand.grid(b_sma = b_sma_vals, r_x1y = r_x1y_vals, r12 = r12_vals)
    parm_combos$b1_ols <- parm_combos$b_sma * parm_combos$r_x1y
  } else {
    parm_combos <- expand.grid(b_sma = b_sma_vals, b1_ols = b1_ols_vals, r12 = r12_vals)
    parm_combos$r_x1y <- parm_combos$b1_ols / parm_combos$b_sma
  }
  
  parm_combos %>% 
    mutate(across(everything(), ~ round(.x, 3))) %>%
    tibble()
}

# NOTE: Cor X1 & Y is positive, Cor X2 & Y is negative, Cor X1 & X2 is negative
parms_mv <- gen_parm_combos_mv(r_x1y_vals = c(.2, .4, .6))
parms_mv <- gen_parm_combos_mv(b1_ols_vals = c(.04, .08, .1))

# TO DO: 
# Potential probs: 1) May be an issue with the 'marginal OLS slope' of Y ~ X1 (vs Y ~ X1 + X2), 2) the residuals from sma_mod may not be appropriate (need to do Peig & Green, 2009 approach?), 

i<- 1
parms_mv[i,]
## sim_SMA_mv() function to simulate the data with specified multivariate relationships
sim_SMA_mv <- function(parms_mv, N = 1000, seed = NULL, beta2 = 0.4) {
  if (!is.null(seed)) set.seed(seed)
  sim_list <- list()
  
  for (i in seq_len(nrow(parms_mv))) {
    cat("Simulating row", i, "\n")
    
    b1_ols <- parms_mv$b1_ols[i]     # desired marginal OLS slope of Y ~ X1
    r12 <- parms_mv$r12[i]          # correlation between X1 and X2
    b_sma <- parms_mv$b_sma[i]     # desired SMA slope
    r_x1y <- parms_mv$r_x1y[i]          # desired correlation between X1 and Y
    
    # Simulate X1 and X2 from bivariate normal with correlation r12
    Sigma <- matrix(c(1, r12, r12, 1), nrow = 2)
    X_vals <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma)
    X1 <- X_vals[, 1]
    X2 <- X_vals[, 2]
    
    # Adjust beta1 to preserve desired marginal OLS slope
    cov12 <- cov(X1, X2)
    var_X1 <- var(X1)
    beta1 <- b1_ols - beta2 * cov12 / var_X1 # cov12 / var_X1 = r12
    
    # Compute residual variance needed to achieve r_x1y and b_ols1
    var_epsilon <- (b1_ols^2 * var_X1) * ((1 - r_x1y^2) / r_x1y^2)
    sd_epsilon <- sqrt(var_epsilon)
    
    # Simulate Y
    eps <- rnorm(N, 0, sd_epsilon)
    Y <- beta1 * X1 + beta2 * X2 + eps
    
    # Fit SMA model and get residuals
    sma_mod <- smatr::sma(Y ~ X1)
    sma_resid <- residuals(sma_mod)
    
    # Store simulation output
    sim_list[[paste0("b_sma=", b_sma, "_r=", r_x1y, "_r12=", r12)]] <- tibble::tibble(
      X1 = X1,
      X2 = X2,
      Y = Y,
      true_r12 = r12,
      true_beta2 = beta2,
      true_b1_ols = b1_ols,
      true_r_x1y = r_x1y,
      true_b1_sma = b_sma,
      est_r12 = cor(X1, X2),
      est_b1_sma = coef(sma_mod)[2],
      est_b2_ols = coef(lm(Y ~ X1 + X2))[3],
      est_b2_noX1 = coef(lm(Y ~ X2))[2],
      est_b2_sma = coef(lm(sma_resid ~ X2))[2]
    )
  }
  dplyr::bind_rows(sim_list, .id = "condition")
}
sim_df_mv <- sim_SMA_mv(parms_mv)

# create_SMA_summary to compare parameter estimates
create_SMA_summary <- function(sim_df_mv) {
  sim_df_mv %>% 
    dplyr::select(starts_with("true"), starts_with("est")) %>%
    distinct()
}
# NOTE:: est_scale_coeff way off, as is the effect of temp from the sma_resid model. 
# Thus far, the est_b2_ols does much better at getting the correct value
create_SMA_summary(sim_df_mv) %>% 
  dplyr::select(-true_beta2)
  #dplyr::select(true_r12, true_b_sma, est_b_sma)

# Plot boxplots to compare simulated parameter value to estimated parm values
plotting_df_mv <- create_SMA_summary(sim_df_mv) %>% 
  pivot_longer(cols = c(est_b2_noX1, est_b2_ols, est_b2_sma), 
               names_to = "Mod_est", values_to = "Temp_b2") %>% 
  mutate(Mod_est = str_remove(Mod_est, "est_b2_"))
plotting_df_mv %>% ggplot(aes(x = Mod_est, y = Temp_b2, color = Mod_est)) + 
  geom_point() +
  geom_boxplot() + 
  geom_hline(yintercept = unique(sim_df_mv$true_beta2), 
             linetype = "dashed", color = "black")

# Do allometric relationships influence the bias in b2 estimates? 
plotting_df_mv %>% ggplot(aes(x = est_b_sma, y = Temp_b2, color = Mod_est)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_hline(yintercept = unique(sim_df_mv$true_beta2), 
             linetype = "dashed", color = "black")

# Depending on this correlation between temp & mass, the pros & cons between multiple regression (as suggested by Ryding, 2022; Freckleton 2002) & allometric residuals (Green, 2001) shift. When there is no correlation between temperature & mass allometric residuals MAY (?) give most unbiased results, but when they are highly correlated it is likely that multiple regression is better? 

# 3 vars manual  --------------------------------------------

# Generate wing length: Strongly decreases with temperature
wing <- 15 - 0.6 * temp + rnorm(n, sd = .5)  

# Generate mass: Decreases with latitude more strongly than wing
mass <- 15 - 5 * temp + .41 * wing + rnorm(n, sd = .5) 

size_temp2 <- data.frame(temp, wing, mass) %>% 
  mutate(wing_log = log(wing), 
         mass_log = log(mass)) %>% 
  tibble()

mod_wm <- sma(wing_log ~ mass_log, data = size_temp2)
mod_mw <- sma(mass_log ~ wing_log, data = size_temp2)
cor(size_temp2$wing_log, size_temp2$mass_log)
summary(mod_wm)
mod_ols <- lm(wing_log ~ mass_log, data = size_temp2)
size_temp_r <- size_temp2 %>% mutate(res_wm = residuals(mod_wm), 
                                     res_mw = residuals(mod_mw),
                                     ols_r = residuals(mod_ols), 
                                     individual = row_number())

## Can flip axes
size_temp_r %>% filter(mass_log > 2.5 & mass_log < 2.6 & res_wm < 0)
size_temp_r %>% filter(individual == 92)

size_temp_r %>% arrange(desc(mass_log)) 
# Wing ~ mass
w_m <- size_temp_r %>%
  ggplot(aes(x = mass_log, y = wing_log)) + 
  geom_point(alpha = .6) +
  geom_point(data = ~filter(.x, individual == 92), 
             size = 5, color = "green") + 
  geom_smooth(method = "lm", color = "red") +
  ggpmisc::stat_ma_line(method = "SMA", color = "blue") 

# Mass ~ wing 
m_w <- size_temp_r %>%
  ggplot(aes(x = wing_log, y = mass_log)) + 
  geom_point(alpha = .6) +
  geom_point(data = ~filter(.x, individual == 92), 
             size = 5, color = "green") + 
  geom_smooth(method = "lm", color = "red") +
  ggpmisc::stat_ma_line(method = "SMA", color = "blue") 

ggarrange(w_m, m_w)

# In SMA -- Residuals are equal & oppositely correlated with X & Y
cor(size_temp_r$mass_log, size_temp_r$res_wm)
cor(size_temp_r$wing_log, size_temp_r$res_wm)

# In OLS -- Residuals are not correlated with X & highly correlated with Y
round(cor(size_temp_r$mass_log, size_temp_r$ols_r), 2) # Not correlated
cor(size_temp_r$wing_log, size_temp_r$ols_r) # Highly correlated

# Comparing residuals
size_temp_r %>% ggplot(aes(x = res_wm, y = ols_r)) + 
  geom_point(alpha = .2) + 
  geom_abline(slope = 1, color = "red") + 
  labs(x = "SMA residuals", y = "OLS residuals")

# Influence of temperature on wingyness
ggplot(data = size_temp_r, aes(x = temp, y = res_wm)) + 
  geom_point(alpha = .6) +
  geom_smooth(method = "lm") + 
  labs(x = "Temperature increase", y = "Wingyness")

# Biologically realistic values -----------------------------------
# Simulate 3 related hypothetical species that have different allometric scaling relationships
# My thought is I can use these hypothetical species to illustrate how size metrics vary depending on the scaling coefficient 
# Is mass:wing ratio the same at different sizes under isometry? 
# What are the implications of taking some sort of residual (e.g., relative wing size, after controlling for allometry) compared to just taking a wing / mass ratio under the different allometric scaling patterns? 

# Sample size & mass values
N <- 1000
Mass <- rnorm(N, mean = 60, sd = 8)
Mass_log <- log(Mass)

# NOTE:: I don't think it is possible to generate data from an exact SMA slope (e.g. 0.33). This is because 1) smatr doesn't include a simulation function, and 2) the correlation is part of what determines b_sma. I.e., b_sma = b_ols / r, so unless we don't include any error term (r = 1) there will always be stochastic variation in the correlation between wing & mass. Thus we must use.. 
## Trial & error:: Identify values of slopes & standard deviations that work to achieve desired b_sma. These values end up in the parms dataframe
# To increase correlation: Increase b_ols or reduce sd of error. Reducing sd of error will create less stochasticity in b_sma 
b_ols <- .18
sd_err <- 0.038 
error <- rnorm(N, 0, sd_err) # Add error on log wing scale
Wing_log <- log(100) + b_ols * Mass_log + error
r <- cor(Wing_log, Mass_log)
r

mod <- summary(lm(Wing_log ~ Mass_log))
b_ols_est <- coef(mod)[2,1]
b_sma <- b_ols / r
b_sma

# Store the combinations of values that worked in a tbl
parms <- tibble(
  Species = c("Isometry", "Hyperallometry", "Hypoallometry"),
  sd = c(.038, .04, .027),
  b_ols = c(.18, .30, .11)
)

## The intercept is critical! Does the intercept matter here? 
log_a <- log(10)
#log_a <- log(1) # Set intercept to 0 

# Use the tbl parms to simulate data
morph_df <- parms %>% rowwise() %>% 
  mutate(Wing_log = list(log_a + (b_ols * Mass_log) + rnorm(N,0,sd))) %>% 
  unnest(Wing_log) %>% 
  mutate(Wing = exp(Wing_log), Mass_log, Mass = exp(Mass_log), 
         r = cor(Mass_log, Wing_log), b_sma = b_ols / r, .by = Species) %>% 
  mutate(Species = factor(Species, levels = c("Hypoallometry", "Isometry", "Hyperallometry"))) %>%
  tibble()
morph_df

## NOTE:: Should really extract the estimated b_ols coefficients to calculate the b_sma values
# lm(wing ~ mass * species, data = morph_df)

# Ensure b_sma comes out as expected
morph_df %>% pull(b_sma) %>% unique()

# Ensure we recover the species specific slopes 
# NOTE:: The elevations are different
mod_sma <- smatr::sma(Wing_log ~ Mass_log * Species, morph_df)
summary(mod_sma)

# Visualize log-log relationship
morph_df %>% 
  ggplot(aes(x = Mass_log, y = Wing_log, color = Species)) + 
  geom_point(alpha = .3) + 
  #geom_smooth(method = "lm", linetype = "dashed", fullrange = TRUE) + 
  ggpmisc::stat_ma_line(method = "SMA", fullrange = TRUE) #+
  #ylim(c(-2, 10)) +
  #xlim(0, 4.6)

# Simple ratios -----------------------------------------------------------
## Understand how simple ratios are affected by scaling theory 
# Unlogged -- Hyperallometric species has the shallowest slope, as expected
morph_df %>% mutate(Mass_wing = Mass / Wing) %>% 
  ggplot(aes(x = Mass, y = Mass_wing, color = Species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

## STILL to do: 
# When log_a = 0, and doing logged mass / wing, you do see flat lines , but I would expect isometry to be flat, hyperallometry to have negative slope, and hypoallometry to have positive slope
morph_df2 <- morph_df %>% 
  mutate(Mass_wing = Mass / Wing, 
         Mass_wing_log = Mass_log / Wing_log)

# NOTE:: According to Jokob (1996), ratios (mass : linear metric) are correlated with body size 
morph_df2 %>% 
  ggplot(aes(x = Mass_log, y = Mass_wing_log, color = Species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", fullrange = TRUE) + 
  labs(x = "Log mass", y = "Mass / wing ratio") 

# Statistical models to compare slopes
mod_sma_mw <- smatr::sma(Mass_wing_log ~ Mass_log * Species, morph_df2)
summary(mod_sma_mw)
# NOTE:: Hypoallometry has steepest slope & hyperallometry has shallowest slope. Makes sense, as mass is the numerator, so in a Hypoallometric species the mass increases relatively more than wing 
summary(lm(Mass_wing_log ~ Mass_log * Species, morph_df2))

# >Same intercepts? --------------------------------------------------------
## How can we obtain the same intercepts, only varying the slopes? 
# Code from GPT, didn't quite work. 

# Define a reference mass for alignment
reference_mass <- mean(Mass_log)

# Use the tbl parms to simulate data
morph_df <- parms %>%
  rowwise() %>%
  mutate(
    Wing_log = list(log_a + (b_ols * Mass_log) + rnorm(N, 0, sd))
  ) %>%
  unnest(Wing_log) %>%
  group_by(Species) %>%
  mutate(
    Mass_log,
    # Adjust wing_log to align intercepts at the reference mass
    Wing_log = Wing_log - (b_ols * reference_mass) + log_a
  ) %>% ungroup()

# Fit the SMA model to check intercepts
mod_sma <- smatr::sma(Wing_log ~ Mass_log * Species, data = morph_df)
summary(mod_sma)

# Visualize to confirm alignment
ggplot(morph_df, aes(x = Mass_log, y = Wing_log, color = Species)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Adjusted Simulation: Same Intercepts Across Species",
       x = "Log(Mass)", y = "Log(Wing)") +
  theme_minimal()

# Understand how simple ratios are affected by scaling theory 
morph_df %>% mutate(Mass_wing = Mass_log / Wing_log) %>% 
  ggplot(aes(x = Mass_log, y = Mass_wing, color = Species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

# Ex1: Temporal bergs - wing^2 / mass ---------------------------------------
# A species appears to adhere to temporal version of Bergmann's rule as both mass and wing decrease w/ temperature, but wing^2 / mass (SA / V) actually decreases with temp (due to hypoallometry). 
# NOTE:: In paper make note that we would know this will happen due to hypoallometry

# Increases in temperature from 1975 to 2025 at different sampling locations 
n <- 1000
temp <- rnorm(n, mean = 1.2, sd = 0.3)  

# Generate wing length: Strongly decreases with temperature
wing <- 15 - 0.6 * temp + rnorm(n, sd = .5)  

# Generate mass: Decreases with latitude, but not as strongly as wing
mass <- 40 - 0.3 * temp + rnorm(n, sd = .5)  

# Compute SA:V ratio (wing^2 / mass)
SA_V <- (wing^2) / mass

# Create data frame
size_temp <- data.frame(temp, wing, mass, SA_V) %>% tibble()

# Check relationships
cor(size_temp)

## Visualize
# How does SA / V ratio change with body size?
size_temp %>% ggplot(aes(x = mass, y = SA_V)) +
  geom_point() + 
  geom_smooth(method = "lm")

p1 <- ggplot(size_temp, aes(x = temp, y = wing)) +
  geom_point(alpha = .4) + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title ="Wing vs Temp", x = "Temperature increase") 

p2 <- ggplot(size_temp, aes(x = temp, y = mass)) +
  geom_point(alpha = .4) + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Mass vs Temp", x = "Temperature increase")

# We would expect that SA:V would increase as temperature increases, but in this case SA:V decreases
p3 <- ggplot(size_temp, aes(x = temp, y = SA_V)) +
  geom_point(alpha = .4) + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "SA:V Ratio vs Temp", x = "Temperature increase")

grid.arrange(p1, p2, p3, nrow = 2)

# Ex2: Neg covariation ---------------------------------------------------
# Individuals are doing different things then the population.. Example, wing & mass show a positive trend with latitude, but actually negatively covary. A migratory bird might diverge in migration strategy & behavior (time-minimizing vs energy-minimizing), where short-distance migrant individuals are fat & short winged, & long-distance migrant individuals are skinny & long-winged.
set.seed(42)

# Generate latitude values
n <- 200  # Number of samples
latitude <- rnorm(n, mean = 50, sd = 10)  # Latitude centered around 50

# Generate wing length: Positively correlated with latitude + some noise
wing <- 20 + 0.1 * latitude + rnorm(n, sd = 1.5)  

# Generate mass: Positively correlated with latitude, but negatively with wing
mass <- 50 + 0.15 * latitude - 0.8 * wing + rnorm(n, sd = 2)  

m_w <- mass/wing

# Create a data frame
bird_data <- data.frame(latitude, wing, mass, m_w)

# Check correlations
cor(bird_data)

# Visualize relationships
p1 <- ggplot(bird_data, aes(x = latitude, y = wing)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Wing vs Latitude")

p2 <- ggplot(bird_data, aes(x = latitude, y = mass)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Mass vs Latitude")

p3 <- ggplot(bird_data, aes(x = wing, y = mass)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Mass vs Wing")

p4 <- ggplot(bird_data, aes(x = latitude, y = m_w)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("")

grid.arrange(p1, p2, p3, nrow = 1) #p4,

# Use allometry? ----------------------------------------------------------
## When should we use allometric scaling theory in estimating body size??

# Let's say you sample from a population, the allometric scaling slope b is estimated from that population. If you resampled, you would get a different slope. So a bird with the exact same wing & mass could fall above the line sometimes, & below the line other times 
# If you sample enough from within a single population, you can probably assume you're getting a pretty good representation of the population & your allometric scaling slope b is probably pretty accurate 
# On other hand, if you are sampling across a continent, or globally, each sample is just a teeny tiny sample from the continental or global population, & your slope would differ significantly if you took a different sample (e.g., moving locations slightly but maintaining the same latitude). 
# A Bayesian framework lends itself well to capture this variability in the estimated slopes

## So it may make sense to use allometric scaling theory in estimating body size (via comparing individuals to the general population) when.. 
# 1) You have high confidence in your SMA line (high sample sizes, or within a single population)
# 2) You want to remove the effect of allometric scaling & keep each individual's position RELATIVE to its expected (mass or wing) given its (wing or mass), given the scaling observed in your empirical sample. When would you or wouldn't you want to do this? 

# >Single population ---------------------------------------------------------
## Extract data from the a single population of a hypothetical species under isometry
iso_spp <- morph_df %>% filter(Species == "Isometry")

# Ensure we recover the expected slope under isometry (0.33)
iso_sma <- smatr::sma(Wing_log ~ Mass_log, iso_spp)
summary(iso_sma)

# Visualize log-log relationship
iso_spp %>% 
  ggplot(aes(x = Mass_log, y = Wing_log)) + 
  geom_point(alpha = .3) + 
  ggpmisc::stat_ma_line(method = "SMA") 

## If we can only sample 300 individuals, how consistent is the SMA slope? 
sma_lines <- map_dfr(1:50, \(rep){
  samp250 <- iso_spp %>% slice_sample(n = 300)
  samp_sma <- smatr::sma(Wing_log ~ Mass_log, samp250)
  tibble(
    pops = "single", rep, 
    int = coef(samp_sma)[1], slope = coef(samp_sma)[2]
    )
})

# Plot variation in lines 
ggplot(data = iso_spp, aes(x = Mass_log, y = Wing_log)) + 
  geom_point(alpha = .1) + 
  geom_abline(data = sma_lines, 
              aes(intercept = int, slope = slope), 
              alpha = .5) +
  ggpmisc::stat_ma_line(method = "SMA", se = FALSE)

# >Mult populations -------------------------------------------------------
## What if we sample 300 individuals from several different populations? 
N_populations <- 100
pop_id <- 1:N_populations
latitude <- round(runif(N_populations, min = 30, max = 50), 2)
N_birds_pp <- 100 # Number of birds per population
b_lat <- 0.8 # Effect of latitude on mass
error_lat <- rnorm(N_populations, 0, 1) # Error associated with latitude

mean_mass <- 28 + (latitude * b_lat) + error_lat
pop_morph_df <- pmap(list(pop_id, latitude, mean_mass), 
                     \(pop, lat, mm){
  mass <- rnorm(N_birds_pp, mm, sd = 8)
  tibble(lat, pop, mass, mass_log = log(mass)) 
}) %>% list_rbind()

# Ensure that the simulated mass values resemble the mass values simulated above in iso_spp (mean 60, sd 8)
iso_spp %>% mutate(mass = exp(Mass_log)) %>% 
  ggplot(aes(x = mass)) +
  geom_density()
pop_morph_df %>% ggplot(aes(x = mass)) +
  geom_density()

# Visualize mass ~ lat in the overall population
pop_morph_df %>% ggplot(aes(x = lat, y = mass)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

## Simulate wing under isometry
iso_parms <- parms %>% filter(Species == "Isometry")

# Generate distinct allometric scaling slopes for each population. This is key to achieve the desired effect
# NOTE:: The variation in the b_ols_pop slopes is what determines the amount of spread when we draw 50x from the global population
b_ols_pop <- rnorm(N_populations, iso_parms$b_ols, .02)
# Create tbl with key parameters for each population
pop_parms <- tibble(pop_id, latitude, b_ols_pop)
# Generate common error values for all populations
error <- rnorm(N_birds_pp, mean = 0, iso_parms$sd)

# Simulate wing values for each population
# IMPORTANT NOTE:: Common intercept with distinct slopes is problematic for maintaining the b_sma value at ~0.33. Can fix this by shifting the estimated slopes down or up a bit (see object slope_dif), but not sure if this is biasing things some other way. Maintaining a common intercept (log(100)) in the generation of wing_log could be problematic... 
pop_morph_df2 <- map2(pop_id, b_ols_pop, \(pop, b_ols){
  morph_df <- pop_morph_df %>% filter(pop == !!pop) #%>% 
    #mutate(mass_log_centered = mass_log - mean(mass_log))
  morph_df %>% mutate(
    wing_log = log(100) + (b_ols * mass_log) + error,
    wing = exp(wing_log)
  ) %>% relocate(wing, .after = mass)
}) %>% list_rbind()

# Visualize relationship between mass & wing
pop_morph_df2 %>% ggplot(aes(x = mass_log, y = wing_log)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

# NOTE:: The overall slope we are recovering is positively biased 
sma_mult_pop <- smatr::sma(wing_log ~ mass_log, data = pop_morph_df2)
coef(sma_mult_pop)[2]

# Adjust all slope estimates by slope_dif so the distribution is maintained but the slopes are centered on top of one another
b_mult_pop <- coef(sma_mult_pop)[2]
b_sing_pop <- coef(iso_sma)[2]
slope_dif <- b_mult_pop - b_sing_pop

# If our 300 individuals are spread over 100 different populations... 
#NOTE:: n can be 300, or take 3 individuals from each population 
sma_lines_pop <- map_dfr(1:50, \(rep){
  samp_pop <- pop_morph_df2 %>% slice_sample(n = 300) # n = 3, by = pop
  samp_sma_pop <- smatr::sma(wing_log ~ mass_log, samp_pop)
  tibble(
    pops = "multi", rep, 
    int = coef(samp_sma_pop)[1], slope = coef(samp_sma_pop)[2] - slope_dif
    )
})

## Plot
# Join with the sma_lines from the single population
sma_lines_compare <- bind_rows(sma_lines, sma_lines_pop)

# Density plot 
ggplot(data = sma_lines_compare, aes(x = slope, color = pops)) + 
  geom_density()

# Plot 50 lines each -- Doesn't do a great job of illustrating the difference
if(FALSE){
  ggplot(data = iso_spp, aes(x = mass_log, y = wing_log)) + 
    geom_point(alpha = .05) + 
    geom_abline(data = sma_lines_pop, 
                aes(intercept = int, slope = slope), 
                color = "red",
                alpha = .5) + 
    geom_abline(data = sma_lines, 
                aes(intercept = int, slope = slope), 
                color = "blue",
                alpha = .5)
}

# SMI ---------------------------------------------------------------------
# Estimate body condition using SMI (Peig & Green, 2009) 
# NOTE:: 3 step process outlined on pages 1886 & 1887

vignette(package = "smatr") # None 

## Calculate the standardized wing index (swi)
mod_swi <- sma(wing_log ~ mass_log, data = pop_morph_df2)

mod_ols <- lm(wing_log ~ mass_log, data = pop_morph_df2)

# L0 is the average mass, essentially allowing for comparison of wing lengths for a given mass 
L0 <- mean(pop_morph_df2$mass)
# Extract the allometric scaling coefficient
b_swi <- coef(mod_swi)[[2]]

# Formula from Peig & Green (2009)
swi_df <- pop_morph_df2 %>% mutate(
  swi = wing * (L0 / mass) ^ b_swi
) %>% arrange(desc(swi))

## Understand example individuals
L0 # Average mass 
mean(pop_morph_df2$wing) # average wing

swi_df %>% slice_head(n = 3) # Large wings relative to their mass
swi_df %>% slice_tail(n = 3) # Small wings relative to their mass

## Plot
swi_df %>% slice_sample(n = 1000) %>% 
  ggplot(aes(x = mass, y = wing, size = swi)) +
  geom_point(alpha = .3)

# >SWI vs wing:mass -------------------------------------------------------
swi_df2 <- swi_df %>% mutate(wing_mass = wing / mass)

# Plot
swi_df2 %>% 
  #pivot_longer() 
  ggplot(aes(x = lat, y = wing_mass)) + # swi
  geom_point(alpha = .1) + 
  geom_smooth(method = "lm")

# These provide similar estimates, but I think in non-isometric species it would be more variable?? That would be pretty cool & important to show! 
summary(lm(wing_mass ~ lat, data = swi_df2))
summary(lm(swi ~ lat, data = swi_df2))


# Thoughts ----------------------------------------------------------------
## Are there cases where we don't need to incorporate a species' allometry? 
# Under isometery (?) simple ratios would produce very similar results as these indices that account for allometric scaling? Andrew: Thinks this makes sense. Under isometry your shape doesn’t change, so ratios should be the same at different sizes. 
# If there IS a slope that is highly different than isometric , need to think critically about whether you want to do some sort of fancy residuals method, or just take a simple ratio

## Next steps: 
# Create a function to simulate data to facilitate quick & easy switches between hypo & hyper allometry, & isometric
# Create a function in 'Mult populations' section and rerun the function 50x to get 

## Recommendations: 
# ALWAYS examine the allometric scaling relationship -- this will help you interpret your results 
# If your species shows allometric scaling near isometry, taking a raw ratio will be nearly identical to directly incorporating allometric scaling theory in estimates of relative mass or wingyness. Or would wing^2 / mass be better?


# Thoughts 2.0 -------------------------------------------------------------
# After meeting with Asad on April 14th
# My simulation was biased towards the OLS method b/c I generated Wing length from the OLS slopes & X1 + X2 + error.
# New plan: Simulate correlated data consistent with Allen's rule & use change in intercepts as a neutral way to measure the true effect of change in temperature on proportions of birds 
## Alternative approaches: 
# 1) Jenna: Do some sort of posterior predictive check to see which model can generate data that looks similar to the observed (i.e., simulated) dataset 
# 2) Asad: Simulate data & run both models, using these parameter estimates as the 'truth'. Then add noise to the simulation & see which method is more robust to additional noise. Essentially, examine the power of each approach to detect the truth 


# >Asad simulation ---------------------------------------------------------
## Simulate data
# Function to generate a symmetric correlation matrix

## OLD
if(FALSE){
  gen_sym_cov <- function(cov_mass_wing, cov_temp_wing, cov_temp_mass){
    matrix(c(
      1, cov_mass_wing, cov_temp_wing,
      cov_mass_wing, 1,   cov_temp_mass,
      cov_temp_wing, cov_temp_mass, 1
    ), nrow = 3)
  }
}

# Desired correlation matrix

# cor_mat <- gen_sym_cov(cov_mass_wing = .3, cov_temp_wing = -.2, cov_temp_mass= -.5) # Produces important point

#cor_mat <- gen_sym_cov(cov_mass_wing = .6, cov_temp_wing = 0, cov_temp_mass= 0)

# Simulate correlated data
#set.seed(1)
#data <- mvrnorm(n = 10000, mu = rep(0, 3), Sigma = cor_mat, empirical = TRUE)

# Transform marginal distributions to be biologically realistic
#Wing       <- log(rescale(data[,1], to = c(40, 50)))# ~ log-Normal(45, 3)
#Mass       <- log(rescale(data[,2], to = c(10, 16))) # ~ log-Normal(13, 2)
#Temp_inc   <- rescale(data[,3], to = c(0.1, 1.6)) # Uniform-like on [0.1, 1.3]

#morph_temp <- tibble(Wing, Mass, Temp_inc)

## 

####



## Generate 3 vars with specified correlations & b_sma value
gen_3var_cov <- function(b_sma_12, r_12 = .5, r_13 = 0, r_23 = 0) {
  S1 <- abs(b_sma_12) 
  S2 <- 1
  cov_12 <- r_12 * S1 * S2
  cov_13 <- r_13 * S1 * 1
  cov_23 <- r_23 * S2 * 1
  matrix(c(
    S1^2,    cov_12, cov_13,
    cov_12,  S2^2,   cov_23,
    cov_13,  cov_23, 1
  ), nrow = 3)
}

parms_mat <- expand_grid(b_sma_12 = c(.22, .33, .44), r_12 = c(.2, .3, .4), r_13 = c(0, -.1, -.2, -.4), r_23 = c(0, -.1, -.2, -.4))

cov_mat <- gen_3var_cov(b_sma_12 = 0.5, r_12 = 0.3, r_13 = -.2, r_23 = -.5)
data <- mvrnorm(n = 10000, mu = rep(0, 3), Sigma = cov_mat, empirical = TRUE)

morph_temp <- tibble(Wing = data[,1], Mass= data[,2], Temp_inc= data[,3])
imap(morph_temp, \(col, name){
  hist(col, main = name)
})
# Check empirical correlations
round(cor(morph_temp), 2)

morph_temp2 <- morph_temp %>% mutate(
  Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)) %>% arrange(Temp_inc) %>% 
  slice_sample(n = 200, by = Temp_bin) %>% 
  filter(!Temp_bin %in% c(1:3, 12:15))

## How can we compare our 2 models to see which one better recovers truth? 
# Extract residuals from SMA
sma_mod <- sma(Wing ~ Mass, data = morph_temp2)
summary(sma_mod)
resid_sma <- residuals(sma_mod)

# Run models to compare the 2 approaches
resid_app <- lm(resid_sma ~ Temp_inc, data = morph_temp2)
summary(resid_app)
ryding_app <- lm(Wing ~ Mass + Temp_inc, data = morph_temp2)
summary(ryding_app)

# Summarize to see how temperature bin maps onto temperature increase
morph_temp2 %>% summarize(mean_temp_inc = mean(Temp_inc),
                          Sample_size = n(), .by = Temp_bin)

# Visualize data 
morph_temp2 %>% ggplot(aes(x = Mass, y = Wing)) + 
  geom_point() + 
  stat_ma_line() + 
  geom_smooth(linetype = "dashed", method = "lm") # Very little difference? 

## Fit SMA models 
# Keep the slope fixed
mod_temp_bin <- sma(Wing ~ Mass + Temp_bin, data = morph_temp2)
coef(mod_temp_bin)
# Allow the slope to vary with binned temp
mod_temp_bin_int <- sma(Wing ~ Mass * Temp_bin, data = morph_temp2)
summary(mod_temp_bin_int)

# Format 'temperature increase' & 'Temp label' columns
format_sma_parms <- function(sma_mod){
  coef(sma_mod) %>% 
    rownames_to_column("Temp_inc") %>% 
    #mutate(Temp_bin = as.numeric(Temp_bin)) %>%
    mutate(
      Temp_inc = str_pad(Temp_inc, side = "left", width = 2, pad = "0"),
      Temp_inc = str_replace(Temp_inc, pattern = "^([0-9])([0-9])$", replacement = "\\1.\\2"),
      Temp_label = paste0(Temp_inc, "°C"), 
      Temp_inc = as.numeric(Temp_inc)) %>%
    tibble() 
}

mod_parms <- format_sma_parms(mod_temp_bin)
mod_parms_int <- format_sma_parms(mod_temp_bin_int)


# Plot Temp increase vs relative wing size (this plot represents TRUTH)
mod_parms_int %>%
  mutate(Temp_inc = as.numeric(Temp_inc)) %>%
  ggplot(aes(x = Temp_inc, y = elevation)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Temperature increase", y = "Relative wing size")

# Residual approach -- We get the correct inference that wingyness increases with temperature
morph_temp2 %>% ggplot(aes(x = Temp_inc, y = resid_sma)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Use predict to generate a new dataframe holding body size constant
newdat <- tibble(Mass = mean(morph_temp2$Mass), Temp_inc = runif(n = 1000, min(morph_temp2$Temp_inc), max(morph_temp2$Temp_inc))) 
pred <- predict(ryding_app, newdata = newdat, interval = 'confidence')
ndf <- cbind(newdat, pred) %>% rename(Wing_pred = fit) %>% 
  tibble()

# Thought it would make sense to see if we recover the slope but doesn't work with no variation in Mass 
#ndf <- ndf %>% mutate(Mass = Mass + rnorm(nrow(ndf), 0, .01))
#sma(Wing_pred ~ Mass, data = ndf)

# We get the incorrect prediction that wingyness decreases with increases in temperature. 
# TO DO: Confirm with Ryding that this is really what she did! 
ndf %>% ggplot(aes(x = Temp_inc, y = Wing_pred)) + 
  geom_point(data = morph_temp2, aes(x = Temp_inc, y = Wing), alpha = .4) + 
  geom_line(color = "blue") +
  geom_ribbon(data=ndf, aes(x=Temp_inc, ymin=lwr, ymax=upr), 
              alpha=0.5, fill = "blue") + 
  labs(y = "Wing_pred (index of wingyness?)")

# EXTRAS ---------------------------------------------------------------------
stop()
# Ex3: Ag vs nat ----------------------------------------------------------
# Can try to simulate this with 1 species, but I think the point will be made more effectively with multiple species. 

# >1 species ---------------------------------------------------------------
set.seed(42)  # For reproducibility

# Generate habitat type (50% AG, 50% NAT)
n <- 200
habitat <- sample(c("AG", "NAT"), size = n, replace = TRUE)

# Baseline values
latitude <- rnorm(n, mean = 50, sd = 10)  # Latitude (not used in model, but realistic)

# Mass: AG birds are fatter
mass <- ifelse(habitat == "AG", 
               40 + rnorm(n, sd = 2),  # AG birds heavier
               36 + rnorm(n, sd = 2))  # NAT birds leaner

# Wing length: AG birds have shorter, fatter wings
wing <- ifelse(habitat == "AG", 
               20 + rnorm(n, sd = 1.5),  # Shorter wings in AG
               23 + rnorm(n, sd = 1.5))  # Longer wings in NAT

# Mass/Wing Ratio: AG birds have a much higher ratio
mass_wing_ratio <- mass / wing

# Create data frame
bird_data <- data.frame(habitat, latitude, mass, wing, mass_wing_ratio)

# Check means per habitat
aggregate(cbind(mass, wing, mass_wing_ratio) ~ habitat, data = bird_data, mean)

# Visualize
p1 <- ggplot(bird_data, aes(x = habitat, y = mass, fill = habitat)) +
  geom_boxplot() + ggtitle("Mass by Habitat")

p2 <- ggplot(bird_data, aes(x = habitat, y = wing, fill = habitat)) +
  geom_boxplot() + ggtitle("Wing Length by Habitat")

p3 <- ggplot(
  bird_data, aes(x = habitat, y = mass_wing_ratio, fill = habitat)
) + geom_boxplot() + 
  ggtitle("Mass/Wing Ratio by Habitat (Stronger Effect)")

grid.arrange(p1, p2, p3, nrow = 2)


# >30 spp ----------------------------------------------------------------
# 1) Ag vs natural habitat: Birds get fatter in agriculture , but also get shorter & fatter wings to increase maneuverability to evade predators. Simulate so you miss key responses when using single variables, but see a marked response when you combine the two metrics. This could replace #2 and just use #2 to cite
# 2) Temporal change through time, missing key responses when using single variables (Weeks, Jirinec, the total body length to mass article)
set.seed(42)  # For reproducibility

# Define parameters
n_species <- 30
n_per_species <- 50  # Individuals per species
habitats <- c("AG", "NAT")

# Generate species names
species_names <- paste0("Species_", 1:n_species)

# Initialize empty data frame
bird_data <- data.frame()
spp_morph <- list()

## Nested loop
# First, loop over each species
for (species in species_names) {
  print(species)
  # Baseline values for wing and mass (species-specific)
  base_mass <- rnorm(1, mean = 35, sd = 7)  # Average mass per species
  base_wing <- rnorm(1, mean = 22, sd = 5)  # Average wing per species
  # Define a 'turn' variable, where 1/2 of the species will not receive the habitat effect for either wing or mass
  turn <- rbinom(2, 1, .5)
  spp_sd_m <- abs(rnorm(1, .5, .5))
  spp_sd_w <- abs(rnorm(1, .3, .2))
  
  # Second, 
  for (i in 1:n_per_species) { 
    print(i)
    habitat <- sample(habitats, 1)  # Randomly assign habitat for each individ
    
    # Habitat effects (individual-specific)
    mass_shift <- rnorm(1, .3, spp_sd_m)
    wing_shift <- rnorm(1, -.8, spp_sd_w) 
    
    # Adjust mass and wing based on habitat
    if (habitat == "AG" & turn[1] == 1) {
      mass <- base_mass + mass_shift + rnorm(1, sd = .1)
    } else if(habitat == "AG" & turn[2] == 1){
      wing <- base_wing + wing_shift + rnorm(1, sd = .1)
    }
    else {
      mass <- base_mass + rnorm(1, sd = spp_sd_m)
      wing <- base_wing + rnorm(1, sd = spp_sd_w)
    }
    
    # Compute mass-to-wing ratio
    mass_wing_ratio <- mass / wing
    
    # Append to data frame
    bird_data <- rbind(bird_data, data.frame(Species = species, Habitat = habitat, Mass = mass, Wing = wing, Mass_Wing_Ratio = mass_wing_ratio)) %>% 
      tibble()
  }
  spp_morph[[species]] <- bird_data
  bird_data <- data.frame() # Reset to an empty dataframe
}

# Scale so coefficients are more comparable
bird_data2 <- spp_morph %>% bind_rows() %>%
  mutate(across(where(is.numeric), scale))
  
# Check the first few rows
head(bird_data2)

# Summary of habitat effects
bird_data2 %>%
  group_by(Habitat) %>%
  summarise(across(c(Mass, Wing, Mass_Wing_Ratio), mean, na.rm = TRUE))

## Bayesian framework 
# NOTE:: This paramterization estimates the mean difference for each habitat value, & then allows each species to have its own intercept & additional variation on the mean habitat. Estimating the overall mean is more computationally efficient than (Habitat | Species)
detach("package:ggpmisc", unload = TRUE)
detach("package:ggpp", unload = TRUE)
DVs <- c("Mass_Wing_Ratio", "Mass", "Wing")
brm_l <- map(DVs, \(DV){
  form <- as.formula(paste0(DV, "~ Habitat + (Habitat | Species)"))
  brms::brm(formula = form,
            family = gaussian(), data = bird_data2,
            chains = 4, iter = 2000, cores = 4)
})
names(brm_l) <- DVs

# Extract fixed effects
fixefs <- map_dfr(brm_l, ~tidy(.x, effects = "fixed"), .id = "DV") %>%
  filter(term == "HabitatNAT") %>%
  select(DV, estimate_fixed = estimate, std.error_fixed = std.error)

# Extract random effects
ranefs <- map_dfr(brm_l, ~tidy(.x, effects = "ran_vals"), .id = "DV") %>%
  filter(term == "HabitatNAT") %>%
  select(DV, Species = level, estimate_ran = estimate, std.error_ran = std.error)

# Join and compute total effect
full_effects <- left_join(ranefs, fixefs) %>%
  mutate(
    Estimate_total = estimate_fixed + estimate_ran,
    Std.error_total = sqrt(std.error_fixed^2 + std.error_ran^2),
    Lower = Estimate_total - (1.96 * Std.error_total), 
    Upper = Estimate_total + (1.96 * Std.error_total), 
  )

## Plot species-specific habitat effects
p1 <- full_effects %>% 
  ggplot(aes(x = Species, y = Estimate_total, color = DV)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  coord_flip() +
  facet_wrap(~DV) + 
  ggtitle("Species-Specific Habitat Effects on Mass/Wing Ratio") +
  ylab("Effect of NAT Habitat on Mass/Wing Ratio") +
  xlab("Species")
p1
# COLOR CODE BASED ON NEGATIVE, POSITIVE, OR OVERLAPS ZERO

## DELETE
# Extract species-specific estimates for each model
species_effects <- map(brm_l, \(mod){
  as.data.frame(brms::ranef(mod)$Species) %>% 
    rownames_to_column("Species")
}) %>% list_rbind(names_to = "DV")

# Clean up data for plotting
species_effects2 <- species_effects %>%
  rename(Estimate = Estimate.HabitatNAT,  # Extract habitat effect
         Lower = Q2.5.HabitatNAT,  # 2.5% credible interval
         Upper = Q97.5.HabitatNAT) %>%  # 97.5% credible interval
  select(Species, DV, Estimate, Lower, Upper)

# >Ex1: Bergs - wing^2 / mass ---------------------------------------
# Original Example 1 using SPATIAL berg's rule 
# A species appears to adhere to spatial version of Bergmann's rule as both mass and wing increase w/ latitude, but wing^2 / mass (SA / V) actually increases with latitude (due to hyperallometry)
n <- 200
latitude <- rnorm(n, mean = 50, sd = 10)  # Latitude centered around 50

# Generate wing length: Strongly increases with latitude
wing <- 15 + 0.2 * latitude + rnorm(n, sd = 1.5)  

# Generate mass: Increases with latitude, but not as strongly as wing
mass <- 40 + 0.1 * latitude + rnorm(n, sd = 2)  

# Compute SA:V ratio (wing^2 / mass)
SA_V <- (wing^2) / mass

# Create data frame
bird_data <- data.frame(latitude, wing, mass, SA_V)

# Check relationships
cor(bird_data)

# Visualize
p1 <- ggplot(bird_data, aes(x = latitude, y = wing)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Wing vs Latitude")

p2 <- ggplot(bird_data, aes(x = latitude, y = mass)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Mass vs Latitude")

p3 <- ggplot(bird_data, aes(x = latitude, y = SA_V)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("SA:V Ratio vs Latitude (Increases)")

grid.arrange(p1, p2, p3, nrow = 2)

## MY ORIGINAL ATTEMPT TO SIMULATE BEFORE CHATGPT
# Sample size & mass values
N <- 1000
mass <- rnorm(N, 60, 8)
sigma <- exp(rnorm(N, 0, .1))
hist(mass)
hist(sigma)

# Define parameters 
a <- 40
scaling_coeff <- setNames(c(.33, .38, .28), 
                          c("Isometry", "Hyperallometry", "Hypoallometry"))


# Generate wing values from mass
morph_df <- imap_dfr(scaling_coeff, \(b, spp){
  wing_log <- log(a) + b * log(mass) + log(sigma)
  wing <- exp(wing_log)
  tibble(species = spp, wing_log, wing, mass, mass_log = log(mass))
})

# Mass & wing on log scale
morph_df %>% ggplot(aes(x = mass_log, y = wing_log, color = species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

# Understand how simple ratios are affected by scaling theory 
morph_df %>% mutate(mass_wing = mass / wing) %>% 
  ggplot(aes(x = mass, y = mass_wing, color = species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")
