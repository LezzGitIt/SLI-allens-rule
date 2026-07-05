## Combined empirical shapeshifting figure
## Layout: Bergmann's (top row) / Inverse Bergmann's (bottom row)
## Mixed-Wingier is excluded from this main figure and reported in the
## Supplementary Information instead (see Scripts/supplementary_info.qmd).
## Reads CSVs produced by Nightjar_shape.R, Weeks_2020_ral.R, Atlantic_birds_shape.R

library(tidyverse)
library(cowplot)
library(ggpubr)
ggplot2::theme_set(theme_cowplot(font_size = 10))

# Load CSVs ---------------------------------------------------------------
parms_all <- bind_rows(
  read_csv("Derived/Csv/Nightjar_parms.csv",  show_col_types = FALSE),
  read_csv("Derived/Csv/Weeks_parms.csv",     show_col_types = FALSE),
  read_csv("Derived/Csv/Atlantic_parms.csv",  show_col_types = FALSE)
) %>%
  filter(
    Direction %in% c("Bergmann's", "Inverse Bergmann's"),
    Approach  %in% c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso"),
    std.error  < 0.5
  ) %>%
  mutate(
    # Replace nightjar common names with scientific names for display
    species   = case_when(
      species == "Nighthawk"      ~ "Chordeiles minor",
      species == "Nightjar"       ~ "Caprimulgus europaeus",
      species == "Whip-poor-will" ~ "Antrostomus vociferus",
      TRUE ~ species
    ),
    Sig_trait = str_to_sentence(Sig_trait),
    Approach  = factor(Approach,
                       levels = c("Ratio", "Ratio2", "Ryding", "Resid_ols", "Sli_est", "Sli_iso"))
  )

# Study aesthetics --------------------------------------------------------
study_order  <- c("Nightjar", "Weeks (2020)", "Atlantic birds")
study_colors <- c("Nightjar" = "#E41A1C", "Weeks (2020)" = "#377EB8", "Atlantic birds" = "#4DAF4A")

approach_labs <- c(
  "Ratio"     = "Wing / Mass",
  "Ratio2"    = "Wing² / Mass",
  "Ryding"    = "Mass as covariate",
  "Resid_ols" = "OLS residuals",
  "Sli_est"   = "SLI estimated",
  "Sli_iso"   = "SLI isometry"
)

shape_scale <- c("Both" = 15, "Mass" = 16, "Neither" = 17, "Wing" = 18)

# Per-direction panel builder ---------------------------------------------
build_direction_plot <- function(df_panel, direction, show_legend = FALSE) {
  species_meta <- df_panel %>%
    group_by(species, Study) %>%
    summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
    mutate(Study = factor(Study, levels = study_order)) %>%
    arrange(Study, mean_est)

  species_levels <- species_meta$species
  axis_colors    <- study_colors[as.character(species_meta$Study)]

  star_df <- df_panel %>%
    group_by(species) %>%
    summarise(y_star          = max(UCI95, na.rm = TRUE),
              rank_consistent = first(rank_consistent),
              .groups = "drop") %>%
    filter(!is.na(rank_consistent) & rank_consistent == TRUE)

  p <- df_panel %>%
    mutate(species = factor(species, levels = species_levels)) %>%
    ggplot(aes(x = species, y = estimate, color = Approach,
               group = interaction(species, Approach))) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = LCI95, ymax = UCI95),
                  alpha = 0.8, width = 0,
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(shape = Sig_trait), size = 1.5,
               position = position_dodge(width = 0.75)) +
    scale_color_hue(labels = approach_labs) +
    scale_shape_manual(values = shape_scale, drop = FALSE) +
    # Approach: 3 cols × 2 rows; Sig_trait: 2 cols × 2 rows — side-by-side via legend.box
    guides(
      color = guide_legend(override.aes = list(size = 3), ncol = 3),
      shape = guide_legend(override.aes = list(size = 3), ncol = 2)
    ) +
    labs(x = NULL, y = expression(beta ~ "on wing shape"), title = direction) +
    theme(
      axis.text.x  = element_text(angle = 60, hjust = 1, size = 7,
                                   color = axis_colors),
      legend.title = element_blank(),
      plot.title   = element_text(size = 10, face = "plain")
    )

  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = species, y = y_star + 0.05, label = "*"),
      color = "black", size = 4, inherit.aes = FALSE
    )
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "top", legend.box = "horizontal")
  }
  p
}

# Build panels (all without legend) ---------------------------------------
directions <- c("Bergmann's", "Inverse Bergmann's")
panels <- map(
  setNames(directions, directions),
  \(d) build_direction_plot(parms_all %>% filter(Direction == d), d)
)

# Shared legend -----------------------------------------------------------
shared_legend <- get_legend(
  build_direction_plot(
    parms_all %>% filter(Direction == "Bergmann's"),
    "Bergmann's",
    show_legend = TRUE
  )
)

# Assemble: legend → Bergmann's → Inverse Bergmann's ----------------------
combined <- plot_grid(
  shared_legend,
  plot_grid(panels[["Bergmann's"]], labels = "a", label_size = 10),
  plot_grid(panels[["Inverse Bergmann's"]], labels = "b", label_size = 10),
  ncol        = 1,
  rel_heights = c(0.08, 1, 1)
)
combined

# Export: full text width × tall enough for three effective rows ----------
n_bergs <- n_distinct(parms_all$species[parms_all$Direction == "Bergmann's"])
fig_width <- max(6, n_bergs * 0.13)   # ~0.13" per species in Bergmann's row

ggsave("Figures/Empirical_combined.png", combined,
       bg = "white", width = fig_width, height = 7.5, units = "in", dpi = 300)

message(sprintf("Saved Figures/Empirical_combined.png  [%.1f\" × 7.5\"]", fig_width))
