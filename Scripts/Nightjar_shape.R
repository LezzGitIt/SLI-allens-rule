## Analysis -- Shapeshifting in nightjars ## 


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(smatr)
library(cowplot)
library(broom)
library(ggpmisc) 
ggplot2::theme_set(theme_cowplot())

load("Rdata/Capri_dfs_07.09.24.Rdata")


# Formatting --------------------------------------------------------------
# Format predictor & responsible variables 
nj_df <- capriA.red2 %>% 
  rename(Wing = Wing.comb, Mass = Mass.comb) %>%
  drop_na(Wing, Mass) %>%
  mutate(B.lat.rd = round(B.Lat, 0),
         log_wing = log(Wing), 
         log_mass = log(Mass),
         wing_mass = Wing / Mass, 
         wing2_mass = Wing^2 / Mass) %>% 
  select(Species, Wing, Mass, B.Lat, log_wing, log_mass, wing_mass, wing2_mass, B.Tavg) 


# Visualize regression approaches -----------------------------------------
# Compare three line-fitting methods
nj_df %>% 
  ggplot(aes(x = log_mass, y = log_wing)) + 
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, color = "red") +
  ggpmisc::stat_ma_line(method = "MA", se = FALSE, color = "orange") +
  ggpmisc::stat_ma_line(method = "SMA", linetype = "dotted", se = FALSE, color = "blue") + 
  facet_wrap(~Species)

# Examine variance in X & Y
nj_df %>% group_by(Species) %>% 
  summarise(var_mass = var(log_mass), var_wing = var(log_wing)) %>%
  mutate(ratio = sqrt(var_wing / var_mass))

# Create list by for each Species
nj_df_l <- nj_df %>% group_split(Species)
names(nj_df_l) <- c("Nighthawk", "Nightjar", "Whip-poor-will")


# SMA models --------------------------------------------------------------
# Run SMA models & extract the residuals 
## NOTE: if you scale first, then the variance of both log_mass & log_wing is 1, & SMA slope = MA slope = 1 × MA slope so these are identical. 
nj_df_l2 <- map(nj_df_l, \(df){
  sma_mod <- sma(log_wing ~ log_mass, data = df, method = "SMA")
  ma_mod <- sma(log_wing ~ log_mass, data = df, method = "MA")
  df %>% mutate(resid_sma = residuals(sma_mod),
                resid_ma = residuals(ma_mod))
})

# Scale by species 
nj_df_l3 <- map(nj_df_l2, \(df){
  df %>% mutate(across(where(is.numeric), scale))
})

# Inspect correlations with body size (mass in this case)
map(nj_df_l3, \(df){
  df %>% summarize(wm_m = cor(wing_mass, Mass),
                   w2m_m = cor(wing2_mass, Mass),
                   resid_m = cor(resid_sma, Mass))
})


# OLS models & extract parms ----------------------------------------------
# Generate Wing mass models, extract parameters via tidy
parms_df <- map(nj_df_l3, \(df){
  mod_coef_sma <- lm(resid_sma ~ B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Resid_sma")
  #mod_coef_ma <- lm(resid_ma ~ B.Tavg, data = df) %>% tidy() %>% 
   # mutate(Approach = "Resid_ma")
  mod_coef_ols <- lm(Wing ~ Mass + B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Ryding")
  mod_coef_ratio <- lm(wing_mass ~ B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass ~ B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Ratio2")
  bind_rows(mod_coef_sma, mod_coef_ols, mod_coef_ratio, mod_coef_ratio2) #, mod_coef_ma
}) %>% list_rbind(names_to = "Species") %>% 
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)


# Plot slope estimates ----------------------------------------------------
# Slope estimates of temperature's impact on shape in nightjars 
legend_labs <- c("Wing / Mass", "Wing² / Mass", "Allometric \nresiduals", "Mass as \ncovariate") #, "Allometric \nresiduals MA"
parms_df %>% filter(term == "B.Tavg") %>% # & !Approach %in% c("Ratio2", "Ryding")) %>% 
  ggplot(aes(x = Species, y = estimate, color = Approach, 
             group = interaction(Species, Approach))) + 
  geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                alpha = .8, width = 0, position = position_dodge(width = 0.75)) +
  geom_point(size = 2, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = NULL, y = "β temperature on body shape") +
  scale_color_hue(labels = legend_labs)
#ggsave("Nightjar_shape.png", bg = "white")
