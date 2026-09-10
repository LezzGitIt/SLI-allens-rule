## Combined empirical shapeshifting figure
## Layout: Bergmann's (top row) / Inverse Bergmann's (bottom row)
## Mixed-Wingier is excluded from this main figure and reported in the
## Supplementary Information instead (see Scripts/supplementary_info.qmd).
## Reads CSVs produced by Nightjar_shape.R, Weeks_2020_ral.R, Atlantic_birds_shape.R

library(tidyverse)
library(cowplot)
library(patchwork)
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

# Significance-classification aesthetics -----------------------------------
# Diverging blue/grey/red so "Sig. longer" and "Sig. stouter" read as opposite
# poles of the same boxplot's zero line, with a neutral grey non-significant
# midpoint (dataviz skill's diverging formula; colors validated with
# scripts/validate_palette.js -- CVD dE 8.7, normal-vision dE 17.8, all >=3:1
# contrast on white). ggplot2's default position_stack() places the *last*
# factor level at the bar's base, so levels are ordered longer -> not sig. ->
# stouter to put "Sig. stouter" at the bottom and "Sig. longer" at the top,
# mirroring the adjoining boxplot's own negative-below/positive-above zero line.
sig_levels <- c("Sig. longer", "Not significant", "Sig. stouter")
sig_colors <- c("Sig. stouter" = "#e34948", "Not significant" = "#898781", "Sig. longer" = "#2a78d6")

# Per-direction panel builder ---------------------------------------------
# Grouped by method (matching fig-compare-approaches / Figure 3's format) rather
# than by species: boxplot summarizes each method's distribution across species,
# jittered points colored by dataset. Percentage in the top-right corner of each
# panel is the share of species whose rank order across methods matches the
# simulation-predicted order (same rank_consistent field used for tbl-species-
# accounting / the Results prose in Allens_methods_sim.qmd's rank-consistency chunk).
build_direction_plot <- function(df_panel, direction) {
  df_panel <- df_panel %>% mutate(Study = factor(Study, levels = study_order))

  pct_consistent <- df_panel %>%
    distinct(species, rank_consistent) %>%
    summarise(pct = round(100 * mean(rank_consistent, na.rm = TRUE))) %>%
    pull(pct)

  df_panel %>%
    ggplot(aes(x = Approach, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Study), width = 0.15, height = 0, alpha = .6, size = 1.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.2,
             fontface = "italic",
             label = paste0(pct_consistent, "% rank-order consistent")) +
    scale_x_discrete(labels = approach_labs) +
    # drop = FALSE so both directions' plots emit an identical 3-entry Study
    # legend (no nightjars are Inverse Bergmann's) and patchwork collects it to one.
    scale_color_manual(values = study_colors, drop = FALSE) +
    # Axis label wording matches fig-compare-approaches (Figure 3, Allens_methods_sim.qmd)
    # so the two figures read as directly comparable quantities.
    labs(x = NULL, y = expression(hat(beta)[T] ~ "on relative appendage length"),
         title = direction, color = NULL) +
    theme(
      axis.text.x  = element_text(angle = 55, hjust = 1, vjust = 1, size = 9,
                                   margin = margin(t = 0)),
      plot.margin  = margin(t = 5.5, r = 5.5, b = 5.5, l = 40),
      legend.position = "top",
      plot.title   = element_text(size = 10, face = "plain")
    )
}

# Per-direction significance-classification panel builder ------------------
# Same species/method denominator as build_direction_plot()'s boxplot (same
# df_panel passed in), so the two panels are directly comparable: for each
# method, the percentage of species whose 95% CI fell entirely above zero
# (Sig. longer), entirely below zero (Sig. stouter), or overlapped zero (Not
# significant).
build_significance_plot <- function(df_panel) {
  df_panel %>%
    mutate(Call = case_when(
      LCI95 > 0 ~ "Sig. longer",
      UCI95 < 0 ~ "Sig. stouter",
      TRUE      ~ "Not significant"
    ),
    Call = factor(Call, levels = sig_levels)) %>%
    count(Approach, Call, .drop = FALSE) %>%
    group_by(Approach) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = Approach, y = pct, fill = Call)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_x_discrete(labels = approach_labs) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = sig_colors, breaks = sig_levels, drop = FALSE) +
    labs(x = NULL, y = "% of species", fill = NULL) +
    theme(
      axis.text.x  = element_text(angle = 55, hjust = 1, vjust = 1, size = 9,
                                   margin = margin(t = 0)),
      plot.margin  = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5),
      legend.position = "top"
    )
}

# Assemble --------------------------------------------------------------
# 2x2 patchwork: boxplot (magnitude/rank order, 2/3 width) + stacked bar (95%-CI
# significance call, 1/3 width) per direction. axes = "collect_x" drops the top
# row's repeated method labels; guides = "collect" pulls the two legends up top.
# patchwork aligns the panel regions, so the Bergmann's (top) and Inverse
# Bergmann's (bottom) data panels come out the same size on their own.
df_berg <- parms_all %>% filter(Direction == "Bergmann's")
df_inv  <- parms_all %>% filter(Direction == "Inverse Bergmann's")

# Only the Bergmann's row carries the legends; the Inverse Bergmann's plots
# suppress theirs so guides = "collect" gathers exactly one of each.
combined <-
  (build_direction_plot(df_berg, "Bergmann's") + labs(tag = "a")) +
  build_significance_plot(df_berg) +
  (build_direction_plot(df_inv, "Inverse Bergmann's") + labs(tag = "b") + guides(color = "none")) +
  (build_significance_plot(df_inv) + guides(fill = "none")) +
  plot_layout(ncol = 2, widths = c(2, 1),
              guides = "collect", axes = "collect_x") &
  theme(legend.position = "top",
        plot.tag = element_text(size = 11, face = "bold"))

combined

fig_width  <- 9.5
fig_height <- 6.8

ggsave("Figures/Empirical_combined.png", combined,
       bg = "white", width = fig_width, height = fig_height, units = "in", dpi = 300)

message(sprintf("Saved Figures/Empirical_combined.png  [%.1f\" × %.1f\"]", fig_width, fig_height))
