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

# Per-direction panel builder ---------------------------------------------
# Grouped by method (matching fig-compare-approaches / Figure 3's format) rather
# than by species: boxplot summarizes each method's distribution across species,
# jittered points colored by dataset. Percentage in the top-right corner of each
# panel is the share of species whose rank order across methods matches the
# simulation-predicted order (same rank_consistent field used for tbl-species-
# accounting / the Results prose in Allens_methods_sim.qmd's rank-consistency chunk).
build_direction_plot <- function(df_panel, direction, show_legend = FALSE) {
  df_panel <- df_panel %>% mutate(Study = factor(Study, levels = study_order))

  pct_consistent <- df_panel %>%
    distinct(species, rank_consistent) %>%
    summarise(pct = round(100 * mean(rank_consistent, na.rm = TRUE))) %>%
    pull(pct)

  p <- df_panel %>%
    ggplot(aes(x = Approach, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Study), width = 0.15, height = 0, alpha = .6, size = 1.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.2,
             fontface = "italic",
             label = paste0(pct_consistent, "% rank-order consistent")) +
    scale_x_discrete(labels = approach_labs) +
    scale_color_manual(values = study_colors) +
    # Axis label wording matches fig-compare-approaches (Figure 3, Allens_methods_sim.qmd)
    # so the two figures read as directly comparable quantities.
    labs(x = NULL, y = expression(hat(beta)[T] ~ "on relative appendage length"), title = direction) +
    theme(
      axis.text.x  = element_text(angle = 55, hjust = 1, vjust = .58, size = 9),
      legend.title = element_blank(),
      plot.title   = element_text(size = 10, face = "plain")
    )

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "top")
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
shared_legend <- ggpubr::get_legend(
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

# Export: fixed size -- x-axis is now the 6 methods (not per-species), so width
# no longer needs to scale with species count the way the old per-species plot did.
fig_width <- 7

ggsave("Figures/Empirical_combined.png", combined,
       bg = "white", width = fig_width, height = 7.5, units = "in", dpi = 300)

message(sprintf("Saved Figures/Empirical_combined.png  [%.1f\" × 7.5\"]", fig_width))
