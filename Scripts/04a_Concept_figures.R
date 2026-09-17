## Two conceptual figures for the manuscript:
##
## Figure 3 (Methods): The Scaled Length Index (SLI), adapted from the Scaled Mass
## Index (Peig & Green 2009).
## (a) General SLI mechanism, illustrated with an independent hypothetical dataset (not
##     one of the paper's simulated species).
## (b) The paper's own two-pathways decomposition, for a simulated hyperallometric example
##     species drawn from the manuscript's simulation grid.
## Saves the combined figure to Figures/SLI_concept.png (embedded in Allens_methods_sim.qmd)
## and the panel-(b) example species to Derived/Rds/Ex_df_hypo_hyper.rds, so fig-ols-sma
## (Discussion) can reuse the same species without regenerating it.
##
## Box 1 (Introduction, "Building intuition for SMA regression"): OLS vs. SMA regression
## under a flip of the X/Y axes, using real Whip-poor-will mass/wing-chord measurements.
## Saves Figures/sma_flip_axes.png (embedded via markdown image in Allens_methods_sim.qmd).

library(tidyverse)
library(cowplot)
library(smatr)
library(grid)
ggplot2::theme_set(theme_cowplot())

source("Scripts/00_Key_allometry_fns.R")
list2env(readRDS("Derived/Rds/simulation_results.rds"), envir = .GlobalEnv)

# Panel (a): general SLI mechanism, illustrated with an independent hypothetical dataset
# (not one of the paper's simulated species -- this panel introduces the standardization
# concept before panel (b) works through the paper's own two-pathways decomposition).
# Exponent and mass spread are exaggerated relative to a typical isometric b=0.33 purely
# for visual clarity of the curvature; SMA is fit to the simulated data rather than using
# the true generating parameters, matching how SLI-estimated works in practice.
set.seed(1)
Sli_concept_n      <- 130
Sli_concept_b_true <- 0.2
Sli_concept_a_true <- 9
Sli_concept_logM   <- rnorm(Sli_concept_n, mean = log(80), sd = 0.5)
Sli_concept_logA   <- log(Sli_concept_a_true) + Sli_concept_b_true * Sli_concept_logM
Sli_concept_df <- tibble(
  M = exp(Sli_concept_logM + rnorm(Sli_concept_n, 0, 0.06)),
  A = exp(Sli_concept_logA + rnorm(Sli_concept_n, 0, 0.06))
)

Sli_concept_fit  <- smatr::sma(log(A) ~ log(M), data = Sli_concept_df)
Sli_concept_bsli <- unname(coef(Sli_concept_fit)[2])
Sli_concept_asli <- exp(unname(coef(Sli_concept_fit)[1]))
Sli_concept_M0   <- mean(Sli_concept_df$M)  # matches sliR::calc_sli(): mean(df[[mass_nm]])

Sli_concept_df <- Sli_concept_df %>%
  mutate(resid = log(A) - (log(Sli_concept_asli) + Sli_concept_bsli * log(M)),
         logM_dist = log(M) - log(Sli_concept_M0))

# Crop the right ~30% of the mass range so the panel isn't dominated by empty space --
# a few of the largest-bodied individuals are simply not shown.
Sli_concept_xmin <- min(Sli_concept_df$M) * 0.9
Sli_concept_xmax_full <- max(Sli_concept_df$M) * 1.06
Sli_concept_xmax <- Sli_concept_xmin + 0.70 * (Sli_concept_xmax_full - Sli_concept_xmin)
Sli_concept_df <- Sli_concept_df %>% filter(M <= Sli_concept_xmax)
Sli_concept_xrng <- c(Sli_concept_xmin, Sli_concept_xmax)

# Two example individuals, balanced: similar |log-mass distance| from M0 (one below, one
# above) and similar |residual| magnitude (one above the population line, one below).
Sli_concept_target_dist <- 0.55 * 0.5 * 2
Sli_concept_small <- Sli_concept_df %>% filter(logM_dist < 0) %>%
  mutate(score = abs(abs(logM_dist) - Sli_concept_target_dist) + abs(abs(resid) - 0.14) * 1.5) %>%
  slice_min(score, n = 1)
Sli_concept_large <- Sli_concept_df %>% filter(logM_dist > 0) %>%
  mutate(score = abs(abs(logM_dist) - Sli_concept_target_dist) + abs(abs(resid) - 0.14) * 1.5) %>%
  slice_min(score, n = 1)
Sli_concept_hl <- bind_rows(Sli_concept_small, Sli_concept_large) %>%
  mutate(id = c("small", "large"),
         a_i = A / M^Sli_concept_bsli,
         A_scaled = a_i * Sli_concept_M0^Sli_concept_bsli)

# Full-width extension for population + both individual curves (dashed), so their shared
# slope is visible throughout, plus the bold arrowed segment marking the actual slide.
Sli_concept_full_curve <- bind_rows(
  tibble(id = "population", a_i = Sli_concept_asli,
         M = seq(Sli_concept_xrng[1], Sli_concept_xrng[2], length.out = 200)),
  Sli_concept_hl %>% reframe(M = seq(Sli_concept_xrng[1], Sli_concept_xrng[2], length.out = 200), .by = c(id, a_i))
) %>%
  mutate(A = a_i * M^Sli_concept_bsli)

Sli_concept_slide <- Sli_concept_hl %>%
  reframe(M = seq(min(M, Sli_concept_M0), max(M, Sli_concept_M0), length.out = 200), .by = c(id, a_i)) %>%
  mutate(A = a_i * M^Sli_concept_bsli)

Sli_concept_arrows <- Sli_concept_hl %>%
  reframe(frac = c(0.3, 0.52, 0.74), .by = c(id, a_i, M)) %>%
  mutate(M_from = M + (Sli_concept_M0 - M) * (frac - 0.08),
         M_to   = M + (Sli_concept_M0 - M) * (frac + 0.08),
         A_from = a_i * M_from^Sli_concept_bsli,
         A_to   = a_i * M_to^Sli_concept_bsli)

Sli_concept_y_top <- max(Sli_concept_hl$A_scaled)
Sli_concept_formula <- "SLI[i] == A[i] ~ bgroup('[', frac(M[0], M[i]), ']')^{b[SLI]}"

p_sli_concept <- ggplot() +
  geom_line(data = Sli_concept_full_curve %>% filter(id == "population"),
            aes(x = M, y = A), color = "grey45", linewidth = 0.7) +
  geom_line(data = Sli_concept_full_curve %>% filter(id != "population"),
            aes(x = M, y = A, group = id), color = "black", linetype = "42", linewidth = 0.5) +
  geom_point(data = Sli_concept_df, aes(x = M, y = A), shape = 1, color = "grey55", size = 1.7, alpha = .85) +
  geom_segment(aes(x = Sli_concept_M0, xend = Sli_concept_M0,
                    y = min(Sli_concept_df$A) * 0.85, yend = Sli_concept_y_top),
               linetype = "dashed", color = "grey20") +
  geom_line(data = Sli_concept_slide, aes(x = M, y = A, group = id), color = "black", linewidth = 1.1) +
  geom_segment(data = Sli_concept_arrows, aes(x = M_from, xend = M_to, y = A_from, yend = A_to),
               arrow = arrow(length = unit(0.22, "cm"), type = "closed"), linewidth = 1.3, color = "black") +
  geom_point(data = Sli_concept_hl, aes(x = M, y = A), size = 3.2, color = "black") +
  geom_point(data = Sli_concept_hl, aes(x = Sli_concept_M0, y = A_scaled), shape = 21, fill = "white",
             color = "black", size = 3.2, stroke = 1.1) +
  annotate("text", x = Sli_concept_M0, y = -Inf, label = "M[0]", parse = TRUE,
           vjust = 1.4, size = 4.2, color = "grey30") +
  annotate("text", x = Sli_concept_xrng[1] + 0.16 * diff(Sli_concept_xrng), y = max(Sli_concept_df$A) * 1.04,
           label = Sli_concept_formula, parse = TRUE, hjust = 0, size = 5.2) +
  scale_x_continuous(limits = Sli_concept_xrng, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = "Mass (M)", y = "Appendage length (A)") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20)),
        plot.margin = margin(t = 12, r = 8, b = 30, l = 5.5))

# Panel (b): the paper's own two-pathways decomposition, for a simulated hyperallometric
# example species (Temp_eff == "Longer") drawn from the manuscript's own simulation grid.
# Exported below to Derived/Rds/Ex_df_hypo_hyper.rds so fig-ols-sma (Discussion) can reuse
# the same species without regenerating it.
Ex_parms_hypo_hyper <- Parms_mat3 %>%
  filter(Scaling != "Isometry" & Temp_eff == "Longer" & r_12 == r_12[1]) %>%
  slice_max(n = 1, order_by = Strength, by = c(Temp_eff, Scaling), with_ties = FALSE)

set.seed(20260710)
Ex_df_hypo_hyper <- gen_ex_data(Ex_parms_hypo_hyper) %>%
  filter(Scaling == "Hyperallometry") %>%
  mutate(Temp_inc = Temp_inc / 3)

saveRDS(Ex_df_hypo_hyper, "Derived/Rds/Ex_df_hypo_hyper.rds")

# SMA and isometric reference lines, computed directly from the slope/intercept
# formulas already introduced in the SMA-intuition callout (Introduction).
r_pathways    <- cor(Ex_df_hypo_hyper$Mass_log, Ex_df_hypo_hyper$Append_log)
sma_slope     <- sign(r_pathways) * sd(Ex_df_hypo_hyper$Append_log) / sd(Ex_df_hypo_hyper$Mass_log)
sma_intercept <- mean(Ex_df_hypo_hyper$Append_log) - sma_slope * mean(Ex_df_hypo_hyper$Mass_log)
iso_slope     <- 1 / 3
iso_intercept <- mean(Ex_df_hypo_hyper$Append_log) - iso_slope * mean(Ex_df_hypo_hyper$Mass_log)

# Expected phenotype for an individual experiencing a representative amount of
# warming (temp_target), used as the shared anchor for both component arrows and the
# marked point below -- rather than an empirical top-quantile average, this ties the
# whole decomposition to one specific, clearly-stated predicted individual.
temp_target <- 0.5
mass_fit  <- lm(Mass_log ~ Temp_inc, data = Ex_df_hypo_hyper)
x_anchor  <- unname(predict(mass_fit, newdata = data.frame(Temp_inc = temp_target)))
y_sma_at_anchor <- sma_intercept + sma_slope * x_anchor
y_iso_at_anchor <- iso_intercept + iso_slope * x_anchor

# Differential-temperature-association displacement: fit each individual's vertical
# distance from the SMA line against its own simulated temperature increase, then
# evaluate that fit at the same temp_target used for x_anchor above, so the red arrow's
# endpoint and the marked point are this same predicted individual's actual position.
Ex_df_hypo_hyper <- Ex_df_hypo_hyper %>%
  mutate(sma_resid = Append_log - (sma_intercept + sma_slope * Mass_log))
resid_fit          <- lm(sma_resid ~ Temp_inc, data = Ex_df_hypo_hyper)
resid_at_target    <- unname(predict(resid_fit, newdata = data.frame(Temp_inc = temp_target)))
y_actual_at_anchor <- y_sma_at_anchor + resid_at_target
predicted_point    <- tibble(Mass_log = x_anchor, Append_log = y_actual_at_anchor)

# The two arrows chain vertically at x_anchor: SMA -> isometric line (green,
# allometry-isometry component), then isometric line -> actual (red,
# temperature-association component), meeting exactly at the isometric line so the
# total span reads as one additive decomposition.
x_line_end  <- max(Ex_df_hypo_hyper$Mass_log)
y_range     <- diff(range(Ex_df_hypo_hyper$Append_log))
x_label     <- min(Ex_df_hypo_hyper$Mass_log)
mean_x      <- mean(Ex_df_hypo_hyper$Mass_log)  # both reference lines cross here by construction

# Green label sits above the isometric line at x_label, in the blank space up there, so its
# connector crosses the isometric line en route to the arrow's midpoint (mirrors the red
# label/connector's use of the blank space above, just closer to the isometric line).
y_green_label <- (iso_intercept + iso_slope * x_label) + 0.14 * y_range
y_red_label   <- y_actual_at_anchor + 0.20 * y_range

# x_blue_label: low on the blue line, in the sparse bottom-right tail of the point cloud
# (off the main correlation diagonal). x_grey_label: where the isometric line has already
# exited the data's own mass range, so text to its right is genuinely blank.
x_blue_label <- mean_x - 0.02
x_grey_label <- x_line_end + 0.02

set.seed(20260818)
Ex_df_pathways_plot <- Ex_df_hypo_hyper %>% slice_sample(n = 700)

p_two_pathways <- ggplot(Ex_df_pathways_plot, aes(x = Mass_log, y = Append_log)) +
  geom_point(alpha = 0.28, size = 1.6, aes(color = Temp_inc)) +
  geom_abline(intercept = sma_intercept, slope = sma_slope, color = "blue", linewidth = 1) +
  geom_abline(intercept = iso_intercept, slope = iso_slope, color = "grey35", linetype = "dashed", linewidth = 1) +
  geom_point(data = predicted_point, aes(x = Mass_log, y = Append_log), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 2.6, stroke = 1) +
  annotate("segment", x = x_anchor, xend = x_anchor, y = y_sma_at_anchor, yend = y_iso_at_anchor,
           color = "forestgreen", linewidth = 1, arrow = arrow(length = unit(0.18, "cm"), ends = "last")) +
  annotate("segment", x = x_anchor, xend = x_anchor, y = y_iso_at_anchor, yend = y_actual_at_anchor,
           color = "firebrick", linewidth = 1, arrow = arrow(length = unit(0.18, "cm"), ends = "last")) +
  annotate("segment", x = x_label, xend = x_anchor, y = y_green_label, yend = (y_sma_at_anchor + y_iso_at_anchor) / 2,
           color = "forestgreen", linewidth = 0.4) +
  annotate("segment", x = x_label, xend = x_anchor, y = y_red_label, yend = (y_iso_at_anchor + y_actual_at_anchor) / 2,
           color = "firebrick", linewidth = 0.4) +
  annotate("text", x = x_label - 0.06, y = y_green_label, label = "Allometry-isometry\ncomponent",
           color = "forestgreen", size = 4.2, hjust = 0.5, fontface = "bold") +
  annotate("text", x = x_label + 0.03, y = y_red_label, label = "Temperature-association\ncomponent",
           color = "firebrick", size = 4.2, hjust = 0.5, fontface = "bold") +
  # "Estimated allometry" is centered (both lines) on the centroid (x=4.4, y=5.11).
  annotate("text", x = 4.45, y = 5.118,
           label = "Estimated allometry", color = "blue", size = 4.2, hjust = 0.5, fontface = "bold") +
  annotate("text", x = 4.45, y = 5.102,
           label = paste0("hat(beta)[SMA] == ", round(sma_slope, 2)), parse = TRUE,
           color = "blue", size = 4.2, hjust = 0.5) +
  # "Isometric reference" sits right of the isometric line, near the top where that line has
  # already exited the point cloud's own mass range -- genuinely blank space, not overlapping.
  annotate("text", x = x_grey_label + 0.12, y = iso_intercept + iso_slope * x_grey_label + 0.018,
           label = "Isometric reference", color = "grey35", size = 4.2, hjust = 0.5, fontface = "bold") +
  annotate("text", x = x_grey_label + 0.12, y = iso_intercept + iso_slope * x_grey_label + 0.004,
           label = "b == 0.33", parse = TRUE, color = "grey35", size = 4.2, hjust = 0.5) +
  scale_color_viridis_c(option = "plasma") +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.26))) +
  labs(x = "Log(mass)", y = "Log(appendage)", color = "Temperature\nincrease") +
  # Legend moves into the blank bottom-right corner of the panel (below "Estimated allometry",
  # right of the point cloud) instead of sitting in its own column outside the panel.
  theme(legend.position = "inside", legend.position.inside = c(0.74, 0.30),
        legend.background = element_blank(), legend.key = element_blank())

# align/axis ensures the two panels' actual plot areas (not just their outer boxes) line
# up in height -- panel (a) has no axis-text chrome and panel (b) does, which otherwise
# leaves their data regions visibly different sizes despite equal cell widths.
Sli_concept_combined <- plot_grid(p_sli_concept, p_two_pathways, labels = c("a", "b"), label_size = 12,
                                   ncol = 2, rel_widths = c(1, 2), align = "h", axis = "tb")
Sli_concept_combined

ggsave("Figures/SLI_concept.png", Sli_concept_combined,
       width = 13, height = 5.7, units = "in", dpi = 300, bg = "white")

# Figure 2 (main text) PDF for sharing outside the pipeline -- see Figs_share/ at repo root.
dir.create("Figs_share", showWarnings = FALSE)
ggsave("Figs_share/Figure2_SLI_concept.pdf", Sli_concept_combined,
       width = 13, height = 5.7, units = "in", bg = "white")

# Box 1 figure (Introduction): OLS vs. SMA regression under a flip of the X/Y axes ------
# Real data (not simulated): Whip-poor-will mass/wing-chord measurements from the same
# Nightjar-family dataset used in 03a_Nightjar_shape.R. Two panels share one data table;
# panel (a) is the natural orientation, panel (b) swaps which variable is on which axis.
# OLS is genuinely refit with the roles swapped in panel (b) -- this is the whole point,
# OLS changes when the axes are exchanged, SMA does not (it is the same line, algebraically
# re-expressed). The same individual is highlighted (in green) in both panels.
Flip_raw <- read_csv("Data/Capri_BA_compare03.29.26.csv", show_col_types = FALSE)
Flip_df <- Flip_raw %>%
  filter(Species == "Whip-poor-will", !is.na(Mass.comb), !is.na(Wing.comb)) %>%
  transmute(Mass = Mass.comb, Append = Wing.comb)

Flip_ols_a <- lm(Append ~ Mass, data = Flip_df)
Flip_a_ols_a <- coef(Flip_ols_a)[1]; Flip_b_ols_a <- coef(Flip_ols_a)[2]
Flip_r <- cor(Flip_df$Mass, Flip_df$Append)
Flip_b_sma_a <- sign(Flip_r) * sd(Flip_df$Append) / sd(Flip_df$Mass)
Flip_a_sma_a <- mean(Flip_df$Append) - Flip_b_sma_a * mean(Flip_df$Mass)

# Highlighted individual: ~average mass (~52g, close to the sample mean) but an unusually
# low wing length -- mirrors the box's original illustrative point (roughly central on the
# X-axis, but far below the cloud in Y), rather than an individual extreme in both variables.
Flip_pt <- Flip_df %>% filter(abs(Mass - 52) < 2) %>% slice_min(Append, n = 1)

Flip_xrng_a <- range(Flip_df$Mass); Flip_yrng_a <- range(Flip_df$Append)
Flip_line_a <- tibble(Mass = seq(Flip_xrng_a[1] - 3, Flip_xrng_a[2] + 3, length.out = 2)) %>%
  mutate(OLS = Flip_a_ols_a + Flip_b_ols_a * Mass, SMA = Flip_a_sma_a + Flip_b_sma_a * Mass)

Flip_panel_a <- ggplot(Flip_df, aes(Mass, Append)) +
  geom_point(color = "grey60", alpha = 0.6, size = 2) +
  geom_line(data = Flip_line_a, aes(y = OLS, color = "OLS"), linewidth = 1) +
  geom_line(data = Flip_line_a, aes(y = SMA, color = "SMA"), linewidth = 1) +
  geom_point(data = Flip_pt, color = "forestgreen", size = 4) +
  scale_color_manual(name = NULL, values = c(OLS = "firebrick", SMA = "steelblue")) +
  labs(x = "Mass (g)", y = "Appendage (mm)") +
  theme_cowplot(font_size = 16) +
  theme(legend.position = "none")

# Panel (b): OLS refit with the roles swapped; SMA algebraically re-expressed (not refit).
Flip_ols_b <- lm(Mass ~ Append, data = Flip_df)
Flip_a_ols_b <- coef(Flip_ols_b)[1]; Flip_b_ols_b <- coef(Flip_ols_b)[2]
Flip_b_sma_b <- 1 / Flip_b_sma_a
Flip_a_sma_b <- -Flip_a_sma_a / Flip_b_sma_a

Flip_yrng_b <- range(Flip_df$Append)
Flip_line_b <- tibble(Append = seq(Flip_yrng_b[1] - 3, Flip_yrng_b[2] + 3, length.out = 2)) %>%
  mutate(OLS = Flip_a_ols_b + Flip_b_ols_b * Append, SMA = Flip_a_sma_b + Flip_b_sma_b * Append)

Flip_panel_b <- ggplot(Flip_df, aes(Append, Mass)) +
  geom_point(color = "grey60", alpha = 0.6, size = 2) +
  geom_line(data = Flip_line_b, aes(y = OLS, color = "OLS"), linewidth = 1) +
  geom_line(data = Flip_line_b, aes(y = SMA, color = "SMA"), linewidth = 1) +
  geom_point(data = Flip_pt, aes(x = Append, y = Mass), color = "forestgreen", size = 4) +
  scale_color_manual(name = NULL, values = c(OLS = "firebrick", SMA = "steelblue")) +
  labs(x = "Appendage (mm)", y = "Mass (g)") +
  theme_cowplot(font_size = 16) +
  theme(legend.position = "none")

# Small legend, drawn manually (cowplot::get_legend() returns an empty grob under the
# currently-installed ggplot2/cowplot combination) and positioned near panel (a)'s own
# label rather than centered/floating at the very top. Sized to match the panels' own
# (enlarged) font size below, not the small default.
Flip_legend <- ggdraw() +
  draw_line(x = c(0.08, 0.15), y = c(0.35, 0.35), color = "firebrick", linewidth = 1.3) +
  draw_label("OLS", x = 0.165, y = 0.35, hjust = 0, size = 16) +
  draw_line(x = c(0.28, 0.35), y = c(0.35, 0.35), color = "steelblue", linewidth = 1.3) +
  draw_label("SMA", x = 0.365, y = 0.35, hjust = 0, size = 16)

# Flip-axes symbol between the two panels: a rounded (curved) double-headed arrow, sized
# and positioned to read as a primary visual element, not a small annotation -- large
# arrow, text about the size of the panels' own axis labels (font_size = 16 below), both
# sitting high in the column, well clear of the "Mass"/"Appendage" axis titles at the
# bottom of panels (a)/(b).
Flip_arrow_grob <- curveGrob(
  x1 = 0.2, y1 = 0.82, x2 = 0.8, y2 = 0.82,
  curvature = -0.5, ncp = 8, square = FALSE,
  arrow = arrow(ends = "both", length = unit(0.16, "inches"), angle = 25),
  gp = gpar(lwd = 3.5)
)
Flip_arrow_panel <- ggdraw() +
  draw_label("flip axes", x = 0.5, y = 0.89, size = 16) +
  draw_grob(Flip_arrow_grob)

# Narrow middle column (just enough room for the arrow) so the two data panels get most
# of the width -- previously a lot of width sat empty between them, leaving each panel (and
# its axis text) smaller than it needed to be.
Flip_top_row <- plot_grid(Flip_panel_a, Flip_arrow_panel, Flip_panel_b, nrow = 1,
                           rel_widths = c(1, 0.2, 1), labels = c("a", "", "b"), label_size = 12)
Flip_combined <- plot_grid(Flip_legend, Flip_top_row, ncol = 1, rel_heights = c(0.06, 1))
Flip_combined

ggsave("Figures/sma_flip_axes.png", Flip_combined,
       width = 10, height = 6.2, units = "in", dpi = 300, bg = "white")

message("Saved Figures/SLI_concept.png, Figures/sma_flip_axes.png, and Derived/Rds/Ex_df_hypo_hyper.rds")
