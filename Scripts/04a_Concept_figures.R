## Figure 3 (Methods): The Standardized Length Index (SLI), adapted from the Scaled Mass
## Index (Peig & Green 2009).
## (a) General SLI mechanism, illustrated with an independent hypothetical dataset (not
##     one of the paper's simulated species).
## (b) The paper's own two-pathways decomposition, for a simulated hyperallometric example
##     species drawn from the manuscript's simulation grid.
## Saves the combined figure to Figures/SLI_concept.png (embedded in Allens_methods_sim.qmd)
## and the panel-(b) example species to Derived/Rds/Ex_df_hypo_hyper.rds, so fig-ols-sma
## (Discussion) can reuse the same species without regenerating it.

library(tidyverse)
library(cowplot)
library(smatr)
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
           label = paste0("b[SMA] == ", round(sma_slope, 2)), parse = TRUE,
           color = "blue", size = 4.2, hjust = 0.5) +
  # "Isometric reference" sits right of the isometric line, near the top where that line has
  # already exited the point cloud's own mass range -- genuinely blank space, not overlapping.
  annotate("text", x = x_grey_label + 0.12, y = iso_intercept + iso_slope * x_grey_label + 0.018,
           label = "Isometric reference", color = "grey35", size = 4.2, hjust = 0.5, fontface = "bold") +
  annotate("text", x = x_grey_label + 0.12, y = iso_intercept + iso_slope * x_grey_label + 0.004,
           label = "b[SMA] == 0.33", parse = TRUE, color = "grey35", size = 4.2, hjust = 0.5) +
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

message("Saved Figures/SLI_concept.png and Derived/Rds/Ex_df_hypo_hyper.rds")
