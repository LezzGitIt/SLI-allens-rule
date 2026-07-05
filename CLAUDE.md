# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an academic R/Quarto project producing a methods paper on Allen's rule and the Standardized Length Index (SLI). The paper argues for using SMA regression over OLS-based approaches when detecting shape-shifting across environmental gradients, and proposes SLI as the preferred metric for measuring relative appendage length. Target journal: Methods in Ecology and Evolution.

**Main manuscript:** `Scripts/qmd/Allens_methods_sim.qmd`
**Supplement:** `Scripts/qmd/supplementary_info.qmd`

Manuscript `.qmd` files live in `Scripts/qmd/`, separate from the analysis `.R` scripts in `Scripts/`. This keeps Quarto's per-render byproducts (`elsarticle.cls`, `elsarticle-harv.bst` — recreated by the `quarto-journals/elsevier` extension's `format-resources` on every render — plus `*_files/`, `*_cache/`) contained to one folder instead of cluttering `Scripts/`. Because `.qmd` files resolve YAML fields (`bibliography:`, `csl:`, `metadata-files:`) and markdown image paths relative to their own location, both files use `../../Suppfiles/...` and `../../Figures/...` (two levels up to project root, not one).

## Rendering Documents

Render a single Quarto document to PDF (execute from project root):
```bash
quarto render Scripts/qmd/<file>.qmd
```

Output PDFs go to `Derived/` (set in `_quarto.yml`), mirroring the source path (e.g. `Derived/Scripts/qmd/<file>.pdf`). The `execute-dir: project` setting means all relative paths in **code chunks** are relative to the project root, not the file's own folder — this is a different resolution rule than the YAML/markdown paths above.

## Key Architecture

### Core functions (`Scripts/Key_allometry_fns.R`)
Must be sourced before running most other scripts. Key functions:
- `gen_data()` / `gen_cov_mat()`: Generate multivariate normal morphological data on the log scale with controllable allometric slope (`b_avg_12`), correlations, and error types (measurement vs. transient)
- `build_sli_slopes_tbl()`: Fit per-group SMA slopes (e.g., by Age × Sex) for SLI estimation; averages slopes across control variables
- `build_group_cor_tbl()`: Per-group mass~appendage OLS correlation table for inspecting allometric relationships
- `calc_sli()`: Compute the Standardized Length Index (Peig & Green 2009); supports per-group SMA slopes via `control` argument
- `calc_lambda()`: Empirical coefficient of variation ratio (var_append / var_mass)
- `classify_direction()`: Classify shapeshifting direction (Bergmann's / Inverse Bergmann's / Mixed / Stable) from tidy lm output per species
- `gen_cor_vars()`: Generate correlated mass/appendage pairs for pairwise exploration
- `format_temp()`: Bin temperature for plotting; `rm_outliers()`: remove >3 SD outliers

### Main manuscript (`Scripts/qmd/Allens_methods_sim.qmd`)
Simulation study comparing six approaches for estimating relative appendage length along a temperature gradient. Approaches: Ratio, Ratio2, Mass-as-covariate (Ryding), OLS residuals, SLI-isometry, SLI-estimated.

### Simulation script (`Scripts/Allometric_scaling_simulation.R`)
Exploratory script for SMA vs OLS slope behaviour under different error structures.

### Empirical scripts
All three export CSVs to `Derived/Csv/` and are combined by `Scripts/Empirical_combined.R`:
- `Scripts/Nightjar_shape.R`: Caprimulgidae museum/banding data; temperature from WorldClim (cached to `Data/Nightjar_temp.rds`); uses `B.Temp` (WorldClim BIO1 at banding location); three species (Nighthawk, Nightjar, Whip-poor-will)
- `Scripts/Weeks_2020_ral.R`: Temporal shape-shifting (1979–2016); year replaces temperature; Wing + Tarsus as appendages; Mass as sole anchor
- `Scripts/Atlantic_birds_shape.R`: Atlantic bird dataset; Wing + Tarsus appendages

### Combined figure (`Scripts/Empirical_combined.R`)
Reads the three CSVs from `Derived/Csv/` and produces `Figures/Empirical_combined.png`. Species colored by Study.

## Key Statistical Concepts

- **SMA vs OLS**: SMA (`smatr` package) assumes error in both X and Y; OLS assumes error only in Y. SMA slope = OLS slope / r, so OLS always underestimates allometric slopes when |r| < 1. For morphometric allometry where mass is measured with error, SMA is more appropriate.
- **SLI (Standardized Length Index)**: `SLI = L × (L0 / M)^b`, where `L0` is the population mean mass and `b` is the SMA slope (Peig & Green 2009). Corrects for body size variation without conditioning on mass in regression.
- **Six methods compared**: (1) Wing/Mass ratio, (2) Wing²/Mass ratio, (3) Mass as covariate (Ryding et al.), (4) OLS residuals of Wing ~ Mass, (5) SLI with isometric slope (b = 0.33), (6) SLI with estimated SMA slope.
- **Shapeshifting direction**: Classified per species as Bergmann's (mass ↓ or wing ↑ with temperature), Inverse Bergmann's, Mixed (both significant), or Stable (neither significant).

## Primary R Packages
`tidyverse`, `smatr`, `cowplot`, `MASS` (mvrnorm), `broom`, `ggpubr`, `MBESS` (cor2cov), `rlang`, `geodata` + `terra` (temperature extraction, Nightjar only)
