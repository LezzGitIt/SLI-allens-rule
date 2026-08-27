# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an academic R/Quarto project producing a methods paper on Allen's rule and the Standardized Length Index (SLI). As of the `reframe-paper1-estimand` merge (session 29), the paper does **not** argue that SMA regression / SLI is superior to OLS-based approaches (ratio, mass-as-covariate, OLS residuals). It instead asks which *estimand* each body-size-standardization method targets, showing that different methods (ratios, OLS-based approaches, SLI-isometry, SLI-estimated) define relative appendage size against different reference relationships and can therefore estimate different quantities — a modelling decision, not neutral preprocessing. The isometric benchmark used to evaluate methods in simulation is presented as one theoretically-motivated reference point, not a universal ground truth every method should be ranked against. Target journal: Methods in Ecology and Evolution.

**Main manuscript:** `Scripts/qmd/Allens_methods_sim.qmd`
**Supplement:** `Scripts/qmd/supplementary_info.qmd`

Because `.qmd` files resolve YAML fields (`bibliography:`, `csl:`, `metadata-files:`) and markdown image paths relative to their own location, both files use `../../Suppfiles/...` and `../../Figures/...` (two levels up to project root, not one).

## Rendering Documents

Render a single Quarto document to PDF (execute from project root):
```bash
quarto render Scripts/qmd/<file>.qmd
```

Quarto's `output-dir: Rendered/` (set in `_quarto.yml`) still mirrors the source path during rendering (e.g. `Rendered/Scripts/qmd/<file>.pdf`) — that's unavoidable, Quarto always mirrors — but a `post-render` hook (`Scripts/flatten_rendered_pdf.R`) copies the finished PDF into a flat `Rendered/<file>.pdf` right after; that's the actual final location to look for output. The nested copy is deliberately *not* deleted at that point — a `pre-render` hook (`Scripts/clean_stale_render_output.R`) removes it instead, right before the *next* render starts. This split exists because the RStudio/Positron Render button serves the just-rendered file from its Quarto-predicted (nested) location via a local preview server immediately after rendering; deleting it in the post-render hook (same render pass) causes a 404 in that preview pane. Deleting it as the first step of the *next* render is safe, since the previous render's preview has already been consumed by then. Both hooks fire on every render (single-file or whole-project, CLI or IDE Render button) via Quarto's project-level hook mechanism, reading `QUARTO_PROJECT_OUTPUT_FILES` to know what to copy. The `execute-dir: project` setting means all relative paths in **code chunks** are relative to the project root, not the file's own folder — this is a different resolution rule than the YAML/markdown paths above.

## Key Architecture

### Core functions (`Scripts/00_Key_allometry_fns.R`)
Must be sourced before running most other scripts. Key functions:
- `gen_data()` / `gen_cov_mat()`: Generate multivariate normal morphological data on the log scale with controllable allometric slope (`b_avg_12`), correlations, and error types (measurement vs. transient)
- `build_sli_slopes_tbl()`: Fit per-group SMA slopes (e.g., by Age × Sex) for SLI estimation; averages slopes across control variables
- `build_group_cor_tbl()`: Per-group mass~appendage OLS correlation table for inspecting allometric relationships
- `calc_sli()`: Compute the Standardized Length Index (Peig & Green 2009); supports per-group SMA slopes via `control` argument
- `calc_lambda()`: Empirical coefficient of variation ratio (var_append / var_mass)
- `classify_direction()`: Classify shapeshifting direction (Bergmann's / Inverse Bergmann's / Mixed / Stable) from tidy lm output per species
- `gen_cor_vars()`: Generate correlated mass/appendage pairs for pairwise exploration
- `format_temp()`: Bin temperature for plotting

### Main manuscript (`Scripts/qmd/Allens_methods_sim.qmd`)
Simulation study comparing six approaches for estimating relative appendage length along a temperature gradient. Approaches: Ratio, Ratio2, Mass-as-covariate (Ryding), OLS residuals, SLI-isometry, SLI-estimated.

### Conceptual figure (`Scripts/04a_SLI_concept_figure.R`)
Builds Figure 3 (main text, Methods): panel (a) the general SLI mechanism on an independent illustrative dataset, panel (b) the two-pathways decomposition on a simulated example species from `simulation_results.rds`. Saves `Figures/SLI_concept.png` (embedded in the manuscript via `png::readPNG()`, matching `04b_Empirical_combined.R`'s pattern) and `Derived/Rds/Ex_df_hypo_hyper.rds` (the panel-(b) example species, reused by `fig-ols-sma` later in the Discussion).

### Simulation script (`Extra_scripts/Allometric_scaling_simulation.R`)
Exploratory script for SMA vs OLS slope behaviour under different error structures. Kept in `Extra_scripts/` alongside the superseded `SMA_body_shape_methods.qmd` (see README.md) since neither is part of the numbered reproducibility pipeline.

### Empirical scripts
All three export CSVs to `Derived/Csv/` and are combined by `Scripts/04b_Empirical_combined.R`:
- `Scripts/03a_Nightjar_shape.R`: Caprimulgidae museum/banding data; temperature from WorldClim (cached to `Data/Nightjar_temp.rds`); uses `B.Temp` (WorldClim BIO1 at banding location); three species (Nighthawk, Nightjar, Whip-poor-will)
- `Scripts/03b_Weeks_2020_ral.R`: Temporal shape-shifting (1979–2016); year replaces temperature; Wing + Tarsus as appendages; Mass as sole anchor
- `Scripts/03c_Atlantic_birds_shape.R`: Atlantic bird dataset; Wing + Tarsus appendages

### Combined figure (`Scripts/04b_Empirical_combined.R`)
Reads the three CSVs from `Derived/Csv/` and produces `Figures/Empirical_combined.png`. Species colored by Study.

## Key Statistical Concepts

- **SMA vs OLS**: SMA (`smatr` package) assumes error in both X and Y; OLS assumes error only in Y. SMA slope = OLS slope / r, so OLS always underestimates allometric slopes when |r| < 1. For morphometric allometry where mass is measured with error, SMA is more appropriate.
- **SLI (Standardized Length Index)**: `SLI = L × (L0 / M)^b`, where `L0` is the population mean mass and `b` is the SMA slope (Peig & Green 2009). Corrects for body size variation without conditioning on mass in regression.
- **Six methods compared**: (1) Wing/Mass ratio, (2) Wing²/Mass ratio, (3) Mass as covariate (Ryding et al.), (4) OLS residuals of Wing ~ Mass, (5) SLI with isometric slope (b = 0.33), (6) SLI with estimated SMA slope.
- **Shapeshifting direction**: Classified per species as Bergmann's (mass ↓ or wing ↑ with temperature), Inverse Bergmann's, Mixed (both significant), or Stable (neither significant).

## Primary R Packages
`tidyverse`, `smatr`, `cowplot`, `MASS` (mvrnorm), `broom`, `ggpubr`, `MBESS` (cor2cov), `rlang`, `geodata` + `terra` (temperature extraction, Nightjar only)
