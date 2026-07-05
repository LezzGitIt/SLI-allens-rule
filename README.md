# Detecting shape-shifting along environmental gradients: a comparison of methods

**Aaron Skinner**

## Overview

This repository contains the analysis code for a methods paper comparing approaches for quantifying relative appendage length (body shape) along environmental gradients, with application to Allen's rule (the ecogeographic pattern of longer appendages in warmer climates).

I compare six commonly used methods — including log-ratio indices, OLS residuals, mass-as-covariate multiple regression, and two Standardized Length Index (SLI) variations — under a simulation framework with controlled allometric structure. In this simulation, the SLI approaches leveraging SMA-estimated slopes is the preferred metric for detecting shape-shifting.

## Repository structure

```         
Scripts/
  qmd/
    Allens_methods_sim.qmd      # Main manuscript (Quarto)
    supplementary_info.qmd      # Supplementary Information (Quarto)
  Key_allometry_fns.R           # Shared functions (sourced by all scripts)
  Run_simulation.R              # Runs the simulation, saves Derived/Rds/simulation_results.rds
  Allometric_scaling_simulation.R  # Exploratory SMA vs OLS simulation
  Nightjar_shape.R              # Empirical case study: Caprimulgidae
  Weeks_2020_ral.R              # Empirical case study: Weeks et al. (2020)
  Atlantic_birds_shape.R        # Empirical case study: Atlantic Forest birds
  Empirical_combined.R          # Combined empirical figure

Suppfiles/                      # Bibliography, journal metadata, title-page partial
_extensions/                    # Quarto elsevier journal-format extension (needed to render)

Data/                           # Raw data (not version-controlled; see below)
Derived/                        # Script outputs (recreatable; not tracked)
Figures/                        # Saved plots (recreatable; not tracked)
```

Manuscript `.qmd` files live in `Scripts/qmd/`, separate from the analysis `.R` scripts, so that Quarto's per-render byproducts don't clutter `Scripts/`. `_quarto.yml`, `_extensions/`, and `Suppfiles/` are tracked because they're required to reproduce the exact PDF output — not just the analysis code.

## Reproducing the analysis

All scripts use paths relative to the project root. Run the simulation and the three empirical case studies before rendering either manuscript:

``` r
source("Scripts/Run_simulation.R")
source("Scripts/Nightjar_shape.R")
source("Scripts/Weeks_2020_ral.R")
source("Scripts/Atlantic_birds_shape.R")
source("Scripts/Empirical_combined.R")
```

Then render the manuscript and supplement from the project root:

``` bash
quarto render Scripts/qmd/Allens_methods_sim.qmd
quarto render Scripts/qmd/supplementary_info.qmd
```

## Data availability

- **Nightjar data**: Data available at: <https://datadryad.org/dataset/doi:10.5061/dryad.pnvx0k6xw>
- **Weeks et al. (2020)**: Weeks, B.C. et al. (2020) Shared morphological consequences of global warming in North American migratory birds. *Ecology Letters* 23, 316–325. Data available at <https://datadryad.org/dataset/doi:10.5061/dryad.8pk0p2nhw>.
- **Atlantic bird traits**: ATLANTIC BIRD TRAITS dataset. Rodrigues et al. (2019) *Ecology* 100, e02647. Data available at **http://onlinelibrary.wiley.com/doi/10.1002/ecy.2647/suppinfo**

## Dependencies

R packages: `tidyverse`, `smatr`, `cowplot`, `patchwork`, `broom`, `MASS`, `MBESS`, `rlang`, `janitor`, `lubridate`, `naniar`

For WorldClim temperature extraction (Nightjar case study only): `geodata`, `terra`

Install all at once:

``` r
install.packages(c("tidyverse", "smatr", "cowplot", "patchwork", "broom",
                   "MASS", "MBESS", "rlang", "janitor", "lubridate", "naniar",
                   "geodata", "terra"))
```

## Citation

\[To be added upon publication\]
