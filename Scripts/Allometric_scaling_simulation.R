# Simulation code to better understand how allometry influences metrics of size & shape shifting
# Goal is to show interesting cases where allometry is important to understanding the true morphological response. Try to... 
# Go beyond just spatial Bergmann's rule (e.g. temporal)
# Go beyond just wing & mass (& even beyond just birds)

## Load libraries
library(tidyverse)
library(janitor)
library(ggplot2)
library(gridExtra)
library(smatr)

# b_sma = b_ols / r -------------------------------------------------------
# Simulate 3 related hypothetical species that have different allometric scaling relationships
# My thought is I can use these hypothetical species to illustrate how size metrics vary depending on the scaling coefficient 
# Is mass:wing ratio the same at different sizes under isometry? 
# What are the implications of taking some sort of residual (e.g., relative wing size, after controlling for allometry) compared to just taking a wing / mass ratio under the different allometric scaling patterns? 

# Sample size & mass values
N <- 1000
mass <- rnorm(N, mean = 60, sd = 8)
mass_log <- log(mass)

# NOTE:: I don't think it is possible to generate data from an exact SMA slope (e.g. 0.33). This is because 1) smatr doesn't include a simulation function, and 2) the correlation is part of what determines b_sma. I.e., b_sma = b_ols / r, so unless we don't include any error term (r = 1) there will always be stochastic variation in the correlation between wing & mass. Thus we must use.. 
## Trial & error:: Identify values of slopes & standard deviations that work to achieve desired b_sma. These values end up in the parms dataframe
# To increase correlation: Increase b_ols or reduce sd of error. Reducing sd of error will create less stochasticity in b_sma 
b_ols <- .18
sd_err <- 0.038 
error <- rnorm(N, 0, sd_err) # Add error on log wing scale
wing_log <- log(100) + b_ols * mass_log + error
r <- cor(wing_log, mass_log)
r

mod <- summary(lm(wing_log ~ mass_log))
b_ols_est <- coef(mod)[2,1]
b_sma <- b_ols / r
b_sma

# Store the combinations of values that worked in a tbl
parms <- tibble(
  species = c("Isometry", "hyperallometry", "hypoallometry"),
  sd = c(.038, .04, .027),
  b_ols = c(.18, .30, .11)
)

# Does the intercept matter here? 
log_a <- log(100)
#log_a <- log(1) # Set intercept to 0 

# Use the tbl parms to simulate data
morph_df <- parms %>% rowwise() %>% 
  mutate(wing_log = list(log_a + (b_ols * mass_log) + rnorm(N,0,sd))) %>% 
  unnest(wing_log) %>% 
  mutate(wing = exp(wing_log), mass_log, mass = exp(mass_log), 
         r = cor(mass_log, wing_log), b_sma = b_ols / r, .by = species)
  
head(morph_df)

## NOTE:: Should really extract the estimated b_ols coefficients to calculate the b_sma values
# lm(wing ~ mass * species, data = morph_df)

# Ensure b_sma comes out as expected
morph_df %>% pull(b_sma) %>% unique()

# Ensure we recover the species specific slopes 
# NOTE:: The elevations are different
mod_sma <- smatr::sma(wing_log ~ mass_log * species, morph_df)
summary(mod_sma)

# Visualize log-log relationship
morph_df %>% 
  ggplot(aes(x = mass_log, y = wing_log, color = species)) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

# Simple ratios -----------------------------------------------------------
## Understand how simple ratios are affected by scaling theory 
# Unlogged -- Hyperallometric species has the shallowest slope, as expected
morph_df %>% mutate(mass_wing = mass / wing) %>% 
  ggplot(aes(x = mass, y = mass_wing, color = species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

## STILL to do: 
# When log_a = 0, and doing logged mass / wing, you do see flat lines , but I would expect isometry to be flat, hyperallometry to have negative slope, and hypoallometry to have positive slope
morph_df2 <- morph_df %>% 
  mutate(mass_wing = mass / wing, 
         mass_wing_log = mass_log / wing_log)

mod_sma_mw <- smatr::sma(mass_wing_log ~ mass_log * species, morph_df2)
summary(mod_sma_mw)

morph_df %>% mutate(mass_wing = mass_log / wing_log) %>% 
  ggplot(aes(x = mass, y = mass_wing, color = species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

# >Same intercepts? --------------------------------------------------------
## How can we obtain the same intercepts, only varying the slopes? 
# Code from GPT, didn't quite work. 

# Define a reference mass for alignment
reference_mass <- mean(mass_log)

# Use the tbl parms to simulate data
morph_df <- parms %>%
  rowwise() %>%
  mutate(
    wing_log = list(log_a + (b_ols * mass_log) + rnorm(N, 0, sd))
  ) %>%
  unnest(wing_log) %>%
  group_by(species) %>%
  mutate(
    mass_log,
    # Adjust wing_log to align intercepts at the reference mass
    wing_log = wing_log - (b_ols * reference_mass) + log_a
  ) %>%
  ungroup()

# Fit the SMA model to check intercepts
mod_sma <- smatr::sma(wing_log ~ mass_log * species, data = morph_df)
summary(mod_sma)

# Visualize to confirm alignment
ggplot(morph_df, aes(x = mass_log, y = wing_log, color = species)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Adjusted Simulation: Same Intercepts Across Species",
       x = "Log(Mass)", y = "Log(Wing)") +
  theme_minimal()

# Understand how simple ratios are affected by scaling theory 
morph_df %>% mutate(mass_wing = mass_log / wing_log) %>% 
  ggplot(aes(x = mass_log, y = mass_wing, color = species)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

# Ex1: Temporal bergs - wing^2 / mass ---------------------------------------
# A species appears to adhere to temporal version of Bergmann's rule as both mass and wing decrease w/ temperature, but wing^2 / mass (SA / V) actually decreases with temp (due to hypoallometry). 
# NOTE:: In paper make note that we would know this will happen due to hypoallometry

# Increases in temperature from 1975 to 2025 at different sampling locations 
n <- 200
temp <- rnorm(n, mean = 1.5, sd = 0.3)  

# Generate wing length: Strongly decreases with temperature
wing <- 15 - 0.6 * temp + rnorm(n, sd = .5)  

# Generate mass: Decreases with latitude, but not as strongly as wing
mass <- 40 - 0.3 * temp + rnorm(n, sd = .5)  

# Compute SA:V ratio (wing^2 / mass)
SA_V <- (wing^2) / mass

# Create data frame
bird_data <- data.frame(temp, wing, mass, SA_V)

# Check relationships
cor(bird_data)

# Visualize
p1 <- ggplot(bird_data, aes(x = temp, y = wing)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title ="Wing vs Temp", x = "Temperature increase") 

p2 <- ggplot(bird_data, aes(x = temp, y = mass)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Mass vs Temp", x = "Temperature increase")

# We would expect that SA:V would increase as temperature increases, but in this case SA:V decreases
p3 <- ggplot(bird_data, aes(x = temp, y = SA_V)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "SA:V Ratio vs Temp (Decreases)", x = "Temperature increase")

grid.arrange(p1, p2, p3, nrow = 2)


# Ex2: Neg covariation ---------------------------------------------------
# Individuals are doing different things then the population.. Example, wing & mass show a positive trend with latitude, but actually negatively covary. A migratory bird might diverge in migration strategy & behavior (time-minimizing vs energy-minimizing), where some individuals are fat & short winged, & others are skinny & long-winged.
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
  ggtitle("Mass vs Wing (Negative Correlation)")

p4 <- ggplot(bird_data, aes(x = latitude, y = m_w)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("")

grid.arrange(p1, p2, p3,  nrow = 2) #p4,



# Use allometry? ----------------------------------------------------------
## When should we use allometric scaling theory in estimating body size??

# Let's say you sample from a population, the allometric scaling slope b is estimated from that population. If you resampled, you would get a different slope. So a bird with the exact same wing & mass could fall above the line sometimes, & below the line other times 
# If you sample enough from within a single population, you can probably assume you're getting a pretty good representation of the population & your allometric scaling slope b is probably pretty accurate 
# On other hand, if you are sampling across a continent, or globally, each sample is just a teeny tiny sample from the continental or global population, & your slope would differ significantly if you took a different sample (e.g., moving locations slightly but maintaining the same latitude). 
# A Bayesian framework lends itself well to capture this variability in the estimated slopes

## So it may make sense to use allometric scaling theory in estimating body size (via comparing individuals to the general population) when.. 
# 1) You have high confidence in your SMA line (high sample sizes, or within a single population)
# 2) You want to remove the effect of allometric scaling & keep each individual's position RELATIVE to its expected (mass or wing) for its (wing or mass), given the scaling observed in your empirical sample. When would you or wouldn't you want to do this? 

# >Single population ---------------------------------------------------------
## Extract data from the a single population of a hypothetical species under isometry
iso_spp <- morph_df %>% filter(species == "Isometry")

# Ensure we recover the expected slope under isometry (0.33)
iso_sma <- smatr::sma(wing_log ~ mass_log, iso_spp)
summary(iso_sma)

# Visualize log-log relationship
iso_spp %>% 
  ggplot(aes(x = mass_log, y = wing_log)) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

## If we can only sample 300 individuals, how consistent is the SMA slope? 
sma_lines <- map_dfr(1:50, \(rep){
  samp250 <- iso_spp %>% slice_sample(n = 50)
  samp_sma <- smatr::sma(wing_log ~ mass_log, samp250)
  tibble(
    pops = "single", rep, 
    int = coef(samp_sma)[1], slope = coef(samp_sma)[2]
    )
})

# Plot variation in lines 
ggplot(data = iso_spp, aes(x = mass_log, y = wing_log)) + 
  geom_point(alpha = .1) + 
  geom_abline(data = sma_lines, 
              aes(intercept = int, slope = slope), 
              alpha = .5)

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
iso_spp %>% mutate(mass = exp(mass_log)) %>% 
  ggplot(aes(x = mass)) +
  geom_density()
pop_morph_df %>% ggplot(aes(x = mass)) +
  geom_density()

# Visualize mass ~ lat in the overall population
pop_morph_df %>% ggplot(aes(x = lat, y = mass)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

## Simulate wing under isometry
iso_parms <- parms %>% filter(species == "Isometry")

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
  samp_pop <- pop_morph_df2 %>% slice_sample(n = 50) # n = 3, by = pop
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
# If your species shows allometric scaling near isometry, taking a raw ratio will be nearly identical to directly incorporating allometric scaling theory in estimates of relative mass or wingyness 

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
  bird_data,  aes(x = habitat, y = mass_wing_ratio, fill = habitat)
) + geom_boxplot() + 
  ggtitle("Mass/Wing Ratio by Habitat (Stronger Effect)")

grid.arrange(p1, p2, p3, nrow = 2)


# >30 spp ----------------------------------------------------------------
# 1) Ag vs natural habitat: Birds get fatter in agriculture , but also get shorter & fatter wings to increase maneuverability to evade predators. Simulate so you miss key responses when using single variables, but see a marked response when you combine the two metrics. This could replace #2 and just use #2 to cite
# 2) Temporal change through time, missing key responses when using single variables (Weeks, Jirinec, the total body length to article)

set.seed(42)  # For reproducibility

# Define parameters
n_species <- 30
n_per_species <- 50  # Individuals per species
habitats <- c("AG", "NAT")

# Generate species names
species_names <- paste0("Species_", 1:n_species)

# Initialize empty data frame
bird_data <- data.frame()

# SIMULATED POORLY - START FROM SCRATCH
# Loop over each species
for (species in species_names) {
  for (i in 1:n_per_species) {
    habitat <- sample(habitats, 1)  # Randomly assign habitat
    
    # Baseline values for wing and mass (species-specific)
    base_mass <- rnorm(1, mean = 35, sd = .5)  # Average mass per species
    base_wing <- rnorm(1, mean = 22, sd = .3)  # Average wing per species
    
    # Habitat effects (species-specific)
    mass_shift <- rnorm(1, 2.5, .1)
    wing_shift <- rnorm(1, -2, .08) 
    
    # Adjust mass and wing based on habitat
    if (habitat == "AG") {
      mass <- base_mass + mass_shift + rnorm(1, sd = .1)
      wing <- base_wing + wing_shift + rnorm(1, sd = .1)
    } else {
      mass <- base_mass + rnorm(1, sd = .1)
      wing <- base_wing + rnorm(1, sd = .1)
    }
    
    # Compute mass-to-wing ratio
    mass_wing_ratio <- mass / wing
    
    # Append to data frame
    bird_data <- rbind(bird_data, data.frame(Species = species, Habitat = habitat, Mass = mass, Wing = wing, Mass_Wing_Ratio = mass_wing_ratio))
  }
}

# Check the first few rows
head(bird_data)

# Summary of habitat effects
bird_data %>%
  group_by(Habitat) %>%
  summarise(across(c(Mass, Wing, Mass_Wing_Ratio), mean, na.rm = TRUE))

## Bayesian framework 
#brm_mw_fe <- brms::brm(Mass_Wing_Ratio ~ Habitat + (1 + Habitat | Species),
#                    family = gaussian(), data = bird_data,
#                    chains = 4, iter = 2000, cores = 4)
#summary(brm_mw_fe)

# Extract species-specific estimates
#species_effects <- as.data.frame(brms::ranef(brm_mw_fe)$Species)

# Clean up data for plotting
species_effects <- species_effects %>%
  mutate(Species = rownames(.),
         Estimate = Estimate.HabitatNAT,  # Extract habitat effect
         Lower = Q2.5.HabitatNAT,  # 2.5% credible interval
         Upper = Q97.5.HabitatNAT)  # 97.5% credible interval

# Plot species-specific habitat effects
ggplot(species_effects, aes(x = Species, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  coord_flip() +
  ggtitle("Species-Specific Habitat Effects on Mass/Wing Ratio") +
  ylab("Effect of NAT Habitat on Mass/Wing Ratio") +
  xlab("Species")


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
