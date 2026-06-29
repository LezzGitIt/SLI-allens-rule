## Examining changes in wing shape in the Brazilian Atlantic Forest ##

# Libraries & data -------------------------------------------------
library(tidyverse)
library(broom)
library(janitor)
library(smatr)
library(cowplot)
library(ggpubr)
library(naniar)
ggplot2::theme_set(theme_cowplot())

Atlantic_birds <- read_csv("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Ecology/Atlantic_bird_traits/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv") 

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

# Inspect -----------------------------------------------------------------
# Possible measurements: Body length, wing, mass, tarsus, tail
Atlantic_birds3 %>% 
  dplyr::select(
   tarsus, wing, tail, mass, body_length, wingspan_mm, head_length_total_mm, bill_width_mm, bill_depth_mm # starts_with("Bill_")
  ) %>%
  naniar::gg_miss_var(show_pct = TRUE)

# 780 species total
Atlantic_birds3 %>% pull(species_) %>% 
  unique()

# Filter dataset --------------------------------------------------------
Atlantic_birds4 <- Atlantic_birds3 %>%
  filter(!is.na(mass) & !is.na(wing) & !is.na(tail)) %>% 
  filter(sex == "Male" & age == "Adult" & year > 1990) %>%
  dplyr::select(species_, wing, mass, tail, log_wing, log_mass, wing_mass, wing2_mass, B.Tavg, latitude, locality, year) 
table(Atlantic_birds4$year)

# Filter species based on latitudinal range and number of observations
Spp_summary <- Atlantic_birds4 %>%
  summarize(min_lat = min(latitude, na.rm = T),
            max_lat = max(latitude, na.rm = T),
            range_lat = min_lat - max_lat,
            num_obs = n(),
            num_locality = length(unique(locality)),
            .by = species_)
Spp_include <- Spp_summary %>% 
  filter(num_obs > 50) %>% 
  arrange(desc(num_obs))

# Create list 
Atl_birds_l <- Atlantic_birds4 %>% 
  right_join(Spp_include) %>% 
  group_split(species_)

# Name list
names(Atl_birds_l) <- map_chr(Atl_birds_l, ~unique(pull(.x, species_)))
length(Atl_birds_l)

# Berg's relationship -------------------------------------------
## To identify possible species to include in shape analysis, first examine the relationship between morphology and temperature. 

# Run model with mass
mass_mod <- map(Atl_birds_l, \(df){
  tidy(lm(mass ~ B.Tavg, data = df)) 
}) %>% list_rbind(names_to = "species_") %>% 
  mutate(dv = "mass") %>% 
  filter(term == "B.Tavg") 

# Run model with wing
wing_mod <- map(Atl_birds_l, \(df){
  tidy(lm(wing ~ B.Tavg, data = df)) 
}) %>% list_rbind(names_to = "species_") %>% 
  mutate(dv = "wing") %>% 
  filter(term == "B.Tavg") 

Bergs <- bind_rows(mass_mod, wing_mod) 

#Bergs %>% filter(species_ == "Conopophaga_melanops")

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
  mutate(sig = p.value < 0.1) %>% 
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

Atlantic_birds4 <- Atlantic_birds3 %>% 
  filter(species_ %in% unique(Spp_keep2$species_))

# Calculate correlations, and SMA and OLS slopes
Slopes_tbl <- Atlantic_birds4 %>% 
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
Atlantic_birds4 %>% 
  ggplot(aes(x = B.Tavg, y = log_mass, color = species_)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  guides(color = "none")

# Compare SMA vs OLS line-fitting methods for allometry 
## NOTE: In general that OLS slopes are near 0
Atlantic_birds4 %>% 
  ggplot(aes(x = log_mass, y = log_wing, color = species_)) + 
  geom_point(alpha = .2) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, alpha = .3) +
  ggpmisc::stat_ma_line(method = "SMA", se = FALSE, alpha = .3) + 
  #facet_wrap(~species_) + 
  guides(color = "none") 

# Allometric correlation --------------------------------------------------
# Some authors note that SMA should only be used if there is a correlation between the two morphological variables. This makes sense because SMA does not distinguish between scatter (correlation) and the functional relationship between X and Y. See Smith (2009) paper. 

# Create list 
Atl_birds_l2 <- Atlantic_birds4 %>% group_split(species_)

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

# 17 species have significant correlation between mass and wing!!
Spp_keep2 <- Spp_metadata2 %>% filter(p_mw < .06) %>% 
  pull(species_)
Spp_keep2

# Keep just species with significant correlation
Atl_birds_l2 <- Atl_birds_l2[Spp_keep2]

# Create shape metrics ----------------------------------------------------
Atl_birds_l3 <- map(Atl_birds_l2, \(df){
  # extract residuals
  ols_mod <- lm(log_wing ~ log_mass, data = df)
  sma_mod <- sma(log_wing ~ log_mass, data = df, method = "SMA")
  est_b_sma <- coef(sma_mod)['slope']
  df <- df %>% mutate(resid_ols = residuals(ols_mod), 
                      resid_sma = residuals(sma_mod)) %>%
    # Calculate SLI
    calc_sli(b_sli = 0.33, Append = wing, rename_col = "sli_isometry") %>% 
    calc_sli(b_sli = est_b_sma, Append = wing, rename_col = "sli_estimated")
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
# Generate wing mass models, extract parameters via tidy
parms_df <- map(Atl_birds_l4, \(df){
  mod_resid_ols <- lm(resid_ols ~ B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Resid_ols")
  #mod_coef_ma <- lm(resid_ma ~ B.Tavg, data = df) %>% tidy() %>% 
  # mutate(Approach = "Resid_ma")
  mod_coef_ols <- lm(wing ~ mass + B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Ryding")
  mod_coef_ratio <- lm(wing_mass ~ B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Ratio")
  mod_coef_ratio2 <- lm(wing2_mass ~ B.Tavg, data = df) %>% tidy() %>% 
    mutate(Approach = "Ratio2")
  mod_sli_iso <- lm(sli_isometry ~ B.Tavg, data = df) %>% tidy() %>%
    mutate(Approach = "Sli_iso")
  mod_sli_est <- lm(sli_estimated ~ B.Tavg, data = df) %>% tidy() %>%
    mutate(Approach = "Sli_est")
  bind_rows(mod_sli_iso, mod_sli_est, mod_resid_ols, mod_coef_ols, mod_coef_ratio, mod_coef_ratio2) #, mod_coef_ma
}) %>% list_rbind(names_to = "species_") %>% 
  mutate(LCI95 = estimate - 1.96 * std.error,
         UCI95 = estimate + 1.96 * std.error)

# Plot slope estimates ----------------------------------------------------
## Slope estimates of temperature's impact on shape 
# Custom label
legend_labs <- c("wing / mass", "wing² / mass", "OLS residuals", "mass as \ncovariate", "SLI estimated", "SLI isometry") #, "Allometric \nresiduals MA"

# Prep plot
Num_spp <- length(unique(parms_df$species_))
parms_df_p <- parms_df %>% filter(term == "B.Tavg") %>%
  left_join(Spp_metadata2[,c("species_", "Direction", "Sig_trait")]) %>% 
  mutate(Approach = factor(Approach),
         species = str_replace(species_, "_", "\n")) 

## GPT
shape_scale <- c("Both" = 15, "Mass" = 16, "Neither" = 17)

plot_shape <- function(df, title = NULL, legend = TRUE, drop_y = FALSE){
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
    ) + scale_shape_manual(values = shape_scale) + 
    scale_color_hue(labels = legend_labs) 
  
  if(drop_y){
    p <- p + theme(
      axis.title.y = element_blank()
    )
  }
  if(!legend) p <- p + theme(legend.position = "none")
  return(p)
}

Shape_plots <- imap(Direction_effect, \(direction, name){
  parms_filt <- parms_df_p %>% filter(Direction == direction)
  drop_y <- direction %in% c("Inverse", "Mixed") # or however you want
  plot <- plot_shape(df = parms_filt, title = name, drop_y = drop_y)
  return(plot)
})




# Custom plotting function 
shape_scale <- c("Both" = 15, "Mass" = 16, "Neither" = 17)
plot_shape <- function(df, title = NULL, legend = TRUE){
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
    ) + scale_shape_manual(values = shape_scale) + 
    scale_color_hue(labels = legend_labs) 
  if(direction %in% c("Inverse", "Mixed")) p <- p + labs(y = NULL)
  if(!legend) p <- p + theme(legend.position = "none")
  return(p)
}

## Map through each direction_effect and plot
# Create character vector for map() function
Direction_effect <- unique(Spp_metadata2$Direction)
Direction_effect <- setNames(Direction_effect, Direction_effect)

# map through
Shape_plots <- imap(Direction_effect, \(direction, name){
  parms_filt <- parms_df_p %>% filter(Direction == direction)
  plot <- plot_shape(df = parms_filt, title = name, dire)
  print(unique(parms_filt$species_))
  return(plot)
})

# Create common legend
com.leg <- get_legend(plot_shape(df = parms_df_p, title = "All"))

# Plot
ggarrange(
  Shape_plots[[1]], Shape_plots[[2]], Shape_plots[[3]], Shape_plots[[4]], 
  common.legend = TRUE, legend.grob = com.leg, labels = "auto"
)

# Save
ggsave("Figures/Atlantic_birds_shape.png", bg = "white", height = 7, width = 9)

## Interpretation
Spp_metadata2 %>% filter(species_ %in% Spp_keep2) %>% 
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
Cvs <- Atlantic_birds4 %>% 
  summarize(N = n(), 
            lambda = calc_lambda(x = mass, y = wing),
            .by = species_) %>%
  mutate(lambda = round(lambda, 2)) %>% 
  arrange(lambda) 

stop()

# Interpretation -------------------------------------------------------------

# This doesn't seem to be a great dataset to test this particular question because there are very few species with >50 observations and a significant correlation between wing and mass. Joliceur 1990 argues that the correlation should be greater than 0.6, which is extremely high. Furthermore, the relationships between temperature and size are pretty small, so temperature may have little effect on size and shape.
#NOTE: I tried with tarsus as well and got similar results

# Next steps --------------------------------------------------------------

# Could try this looking at how shape changes with the amount of forest or the habitat where the bird was captured
# Could try with some sort of bill measurement