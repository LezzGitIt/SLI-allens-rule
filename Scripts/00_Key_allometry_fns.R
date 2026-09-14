## Relative appendage length key functions script

### The simulation and SLI functions now live in the sliR package. Install with remotes::install_github("LezzGitIt/sliR@v0.1.0").
# Called directly as sliR::... at the analysis-script call sites: calc_sli(), build_sli_slopes_tbl(), sim_allometric() (see gen_ex_data below), sim_correlated().
# Kept as thin local wrappers here, because the paper's vocabulary differs from sliR's generic API:
#   gen_data()   -> sim_allometric()  (adds the paper's Temp_inc/Temp_bin columns via format_temp(); driven by pmap() over a b_avg_12/r_12/r_13/r_23 parameter grid whose column names must match the wrapper's arguments)
#   gen_cov_mat()-> build_cov_mat()   (rescales the gradient block to sd_temp; displayed as teaching content in Extra_scripts/SMA_body_shape_methods.qmd)
#   build_group_cor_tbl() -> its sliR namesake (renames r/p_value back to the r_mw/p_mw that ~6 downstream filters per empirical script depend on)
# Still local, deliberately not in sliR: format_temp, run_sma_mod, format_sma_parms, gen_ex_data, calc_lambda, classify_direction, build_sli_mass_cor_tbl, test_group_effect, calc_sli_hierarchical.

# Load required libraries
# MASS is no longer used by this file, but is left attached because supplementary_info.qmd sources this script without loading MASS itself; dropping it here would change that document's search path.
library(MASS)
library(tidyverse)
library(sliR)

### The wrappers below were validated for statistical equivalence against sliR 0.1.0. Fail loudly if an older sliR is installed, so a stale package cannot silently change manuscript results. Install/update with: remotes::install_github("LezzGitIt/sliR@v0.1.0").
stopifnot(
  "sliR (>= 0.1.0) is required; install with remotes::install_github('LezzGitIt/sliR@v0.1.0')" =
    utils::packageVersion("sliR") >= "0.1.0"
)

# Creation of temperature bins for plotting and examination of scaling intercepts
format_temp <- function(df){
  df %>% mutate(
    Temp_bin = cut(Temp_inc, breaks = 15, labels = FALSE, ordered_result = TRUE)) %>%
    arrange(Temp_inc) %>%
    mutate(Temp_bin = case_when(
      Temp_bin %in% c(1:5) ~ 5,  
      Temp_bin %in% c(11:15) ~ 11,
      .default = Temp_bin
    )) #%>% slice_sample(n = 200, by = Temp_bin) 
}
?cut_number # Consider cut_number to make groups with equal number of individuals

# Spell out a small integer (0-99) as a capitalized English word, for inline reporting
# where a live `r`-computed value happens to start a sentence (numerals shouldn't open
# a sentence in running prose, but the value itself should stay reproducible rather
# than being hardcoded as a literal word).
number_to_word <- function(n) {
  stopifnot(n >= 0, n <= 99, n == round(n))
  ones <- c("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
            "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen",
            "seventeen", "eighteen", "nineteen")
  tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety")
  word <- if (n < 20) {
    ones[n + 1]
  } else {
    remainder <- n %% 10
    ten_word <- tens[(n %/% 10) - 1]
    if (remainder == 0) ten_word else paste0(ten_word, "-", ones[remainder + 1])
  }
  paste0(toupper(substring(word, 1, 1)), substring(word, 2))
}

## Generate data (on the log-scale) using var-cov matrix
# Specify measurement error and transient 'error'
### Now a thin wrapper over sliR::sim_allometric(). Arguments and output columns are unchanged, so every call site still works. Temperature is drawn as a normal gradient, matching the multivariate-normal draw this function used to do; sliR's own default gradient is uniform.
# vary_sd was previously hidden inside gen_cov_mat(), where it silently overwrote whatever sd_log_morph the caller passed. It is surfaced here as an argument, defaulting to TRUE to preserve the historical behaviour. Setting vary_sd = FALSE honours sd_log_morph. It affects raw-scale dispersion and skew only: the log-scale slopes and correlations are invariant to it.
gen_data <- function(n = 3000,
                     b_avg_12 = 0.33,
                     r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                     mean_mass = 80, mean_append = 180, mean_temp = 1,
                     sd_log_morph = 0.10,
                     sd_temp = 0.18,
                     vary_sd = TRUE,
                     meas_error = 0,
                     transient_error_mass = 0,
                     transient_error_append = 0) {

  if (vary_sd) sd_log_morph <- runif(1, 0.05, 0.09)

  sim <- sliR::sim_allometric(
    n             = n,
    b_avg         = b_avg_12,
    r_app_mass    = r_12,
    sd_log_mass   = sd_log_morph,
    mean_append   = mean_append,
    mean_mass     = mean_mass,
    gradient      = "Temp_inc",
    gradient_dist = "normal",
    mean_gradient = mean_temp,
    sd_gradient   = sd_temp,
    r_grad_app    = r_13,
    r_grad_mass   = r_23,
    meas_error             = meas_error,
    transient_error_append = transient_error_append,
    transient_error_mass   = transient_error_mass,
    trim_sd = 3
  )

  format_temp(sim)
}

## Generate the variance-covariance matrix (on the log scale), for log(appendage), log(mass) and temperature
### Now a thin wrapper over sliR::build_cov_mat(). The allometry has two degrees of freedom, so b_avg_12 and r_12 together pin both slopes, and one standard deviation then fixes the scale; sliR solves for the other. Output is an unnamed 3x3 matrix, as before.
# vary: which morphological trait sd_log_morph refers to. The two standard deviations covary, so naming one solves for the other.
# vary_sd: TRUE draws sd_log_morph from Uniform(0.05, 0.09) and ignores whatever was passed. Kept for backwards compatibility. It affects raw-scale dispersion and skew only: the log-scale slopes and correlations are invariant to it.
gen_cov_mat <- function(b_avg_12 = 0.33,
                        r_12 = 0.3, r_13 = -0.1, r_23 = -0.1,
                        vary = c("mass", "append"),
                        sd_log_morph = .07,
                        vary_sd = TRUE,
                        sd_temp = 0.18) {

  vary <- match.arg(vary)
  if (vary_sd) sd_log_morph <- runif(1, 0.05, 0.09)

  sd_arg <- if (vary == "mass") list(sd_log_mass = sd_log_morph) else list(sd_log_append = sd_log_morph)

  Sigma <- do.call(sliR::build_cov_mat, c(
    list(b_avg = b_avg_12, r_app_mass = r_12,
         gradient = "Temp", r_grad_app = r_13, r_grad_mass = r_23),
    sd_arg
  ))

  ## build_cov_mat() always places the gradient at unit SD, because it describes a correlation structure rather than temperature's units. Rescale that block to sd_temp.
  scale_temp <- diag(c(1, 1, sd_temp))
  unname(scale_temp %*% Sigma %*% scale_temp)
}


# Fit an SMA model of Append_log ~ Mass_log, optionally allowing the slope to
# vary with binned temperature (Temp_bin). Used to validate simulated species'
# direction of shape change via their SMA intercepts.
run_sma_mod <- function(df, interaction = FALSE) {
  if (!interaction) {
    smatr::sma(Append_log ~ Mass_log + Temp_bin, data = df, method = "SMA")
  } else {
    smatr::sma(Append_log ~ Mass_log * Temp_bin, data = df, method = "SMA")
  }
}

# Tidy the per-Temp_bin coefficients (intercept/slope) from an sma() model
# fit with run_sma_mod(), decoding the bin label back to a numeric/label pair.
format_sma_parms <- function(sma_mod) {
  coef(sma_mod) %>%
    tibble::rownames_to_column("Temp_inc") %>%
    dplyr::mutate(
      Temp_inc = stringr::str_pad(Temp_inc, side = "left", width = 2, pad = "0"),
      Temp_inc = stringr::str_replace(Temp_inc, "^([0-9])([0-9])$", "\\1.\\2"),
      Temp_label = paste0(Temp_inc, "°C"),
      Temp_inc = as.numeric(Temp_inc)
    ) %>%
    tibble::tibble()
}

# Regenerate raw individual-level data for one or more hypothetical species
# (rows of a Parms_mat-style tibble) via gen_data(), for illustrative figures.
gen_ex_data <- function(Parms_mat, transient_error_mass = 0, transient_error_append = 0) {
  Cols <- Parms_mat %>% dplyr::select(dplyr::starts_with(c("b_", "r_")))
  Parms_mat %>%
    dplyr::mutate(coefs = purrr::pmap(Cols, \(...) gen_data(...,
                                              transient_error_mass   = transient_error_mass,
                                              transient_error_append = transient_error_append))) %>%
    tidyr::unnest(coefs)
}

# Per-group mass ~ appendage OLS correlation table.
# Returns one row per group combination (age × sex) with n, b_ols, r_mw, and p_mw.
# Printed for user inspection: groups with r_mw < cor_min or p_mw >= threshold
# lack a meaningful allometric relationship and should not drive per-group SMA slopes.
### Now a thin wrapper over sliR::build_group_cor_tbl(), which names the correlation and its p-value r and p_value, and returns an ungrouped tibble. Both are translated back here (r_mw, p_mw, grouped by the control variables) so existing call sites and their downstream filters are unaffected.
build_group_cor_tbl <- function(df, Append, Mass = Mass, control,
                                 unknown_codes = c("Unk", "U", "Unknown")) {
  sliR::build_group_cor_tbl(df, Append = {{ Append }}, Mass = {{ Mass }},
                            control = control, unknown_codes = unknown_codes) %>%
    dplyr::rename(r_mw = r, p_mw = p_value) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(control)))
}

# Per-species Pearson correlation between wing-based SLI-isometry/SLI-estimated
# (columns must be named sli_isometry/sli_estimated, as produced by calc_sli())
# and body mass. Exported alongside the ratio-mass correlations for tbl-ratio-mass-summary.
# SLI-estimated contributes no row when it is NA for every individual (species
# failed the upstream per-group allometric-correlation filter).
build_sli_mass_cor_tbl <- function(df, Mass = Mass) {
  mass_vec <- rlang::eval_tidy(rlang::enquo(Mass), df)
  ct_iso <- cor.test(df$sli_isometry, mass_vec)
  rows <- tibble::tibble(Metric = "Sli_iso", n = nrow(df),
                         r = as.numeric(ct_iso$estimate), p_value = ct_iso$p.value)
  if (any(!is.na(df$sli_estimated))) {
    ok <- !is.na(df$sli_estimated)
    ct_est <- cor.test(df$sli_estimated[ok], mass_vec[ok])
    rows <- dplyr::bind_rows(rows, tibble::tibble(Metric = "Sli_est", n = sum(ok),
                                                   r = as.numeric(ct_est$estimate), p_value = ct_est$p.value))
  }
  rows
}

# Test whether a grouping variable (Age/Sex) significantly affects a morphometric
# DV, controlling for the environmental gradient, per species in df_list. Drops
# unknown-coded/NA iv rows, then individual groups below min_n_per_group; returns
# NULL for a species with <2 valid groups remaining (too few to test/use as a
# covariate). Used to decide, per species, whether Age/Sex should be considered
# as covariates at all -- shared across the three empirical scripts so the same
# significance-gating logic can't drift between them.
test_group_effect <- function(df_list, dv, iv, gradient, min_n_per_group = 0,
                              unknown_codes = c("Unk", "U", "Unknown")) {
  purrr::map(df_list, \(df) {
    df_filt <- df %>% dplyr::filter(!is.na(.data[[iv]]) & !(.data[[iv]] %in% unknown_codes))
    if (nrow(df_filt) < 10) return(NULL)

    valid_grps <- df_filt %>%
      dplyr::count(.data[[iv]]) %>%
      dplyr::filter(n >= min_n_per_group) %>%
      dplyr::pull(.data[[iv]])
    if (length(valid_grps) < 2) return(NULL)
    df_filt <- df_filt %>% dplyr::filter(.data[[iv]] %in% valid_grps)

    n_counts <- df_filt %>%
      dplyr::count(.data[[iv]], name = "n") %>%
      dplyr::mutate(lbl = paste0("n_", tolower(.data[[iv]]))) %>%
      dplyr::select(lbl, n) %>%
      tidyr::pivot_wider(names_from = lbl, values_from = n)
    fmla <- as.formula(paste(dv, "~", gradient, "+", iv))
    broom::tidy(lm(fmla, data = df_filt)) %>%
      dplyr::filter(stringr::str_starts(term, iv)) %>%
      dplyr::mutate(dv = dv) %>%
      dplyr::bind_cols(n_counts)
  }) %>%
    purrr::list_rbind(names_to = "species_")
}

# Fits SLI-estimated's exponent (b_sli) per individual via a reliability-gated cascade,
# finest to coarsest:
#   (1) the individual's own combined cell across all `covariates` (only possible/attempted
#       when 2 covariates are supplied -- e.g. a species' own Age x Sex combination),
#   (2) the individual's own single-covariate marginal group(s) (own Age class and/or own
#       Sex class, pooling across the other covariate) -- averaged together when both pass
#       reliability, used as-is (unblended) when only one does.
#   (3) the species-wide pooled slope (fit on every individual passed in, regardless of
#       covariate class) -- reached when neither marginal passes (or, for a species with
#       only one significant covariate, when its own marginal fails), and also the value
#       used for every individual when `covariates` is empty.
# A level is used only when it clears BOTH a minimum sample size and a reliable mass~
# appendage correlation (r >= cor_min, p < cor_p_max); otherwise the cascade drops to the
# next-coarsest level. Individuals with unrecorded Age/Sex are collapsed to their own "Unk"
# class and cascade through the identical levels -- not a special case or an automatic
# fallback to pooled.
# Species-level gate is internal, not the caller's responsibility: the cascade is only
# attempted at all if the species-wide (all-individuals) mass~appendage correlation is
# itself reliable; otherwise sli_estimated = NA for every individual, even ones that would
# have had a reliable cell or marginal group -- this keeps SLI-estimated defined for either
# all or none of a species' individuals, matching the other five approaches' shared N.
calc_sli_hierarchical <- function(df, Append, Mass = Mass, covariates = character(0),
                                   n_min_cell = 50, n_min_marginal = 100,
                                   cor_min = 0.3, cor_p_max = 0.05,
                                   unknown_codes = c("Unk", "U", "Unknown"),
                                   M0 = NULL, rename_col = "sli_estimated") {
  app_q   <- rlang::enquo(Append)
  mass_q  <- rlang::enquo(Mass)
  app_nm  <- rlang::as_label(app_q)
  mass_nm <- rlang::as_label(mass_q)
  stopifnot(length(covariates) %in% 0:2)

  if (is.null(M0)) M0 <- mean(df[[mass_nm]], na.rm = TRUE)

  d <- df %>%
    dplyr::mutate(.log_app = log(.data[[app_nm]]), .log_mass = log(.data[[mass_nm]]))

  # Fits an SMA slope (log scale, as allometric slopes require) plus mass~appendage
  # reliability stats (raw scale, matching the cor_min/cor_p_max convention already used
  # everywhere else in this pipeline -- Spp_keep_vec, build_group_cor_tbl's r_mw/p_mw --
  # so the same species/group is never judged reliable by one check and not the other).
  fit_group <- function(sub) {
    if (nrow(sub) < 3) return(tibble::tibble(n = nrow(sub), r = NA_real_, p = NA_real_, slope = NA_real_))
    r <- suppressWarnings(stats::cor(sub[[app_nm]], sub[[mass_nm]]))
    p <- tryCatch(stats::cor.test(sub[[app_nm]], sub[[mass_nm]])$p.value, error = \(e) NA_real_)
    slope <- tryCatch(
      unname(stats::coef(smatr::sma(.log_app ~ .log_mass, data = sub, method = "SMA"))["slope"]),
      error = \(e) NA_real_)
    tibble::tibble(n = nrow(sub), r = r, p = p, slope = slope)
  }
  reliable <- function(tbl, n_min) !is.na(tbl$r) & tbl$n >= n_min & tbl$r >= cor_min & !is.na(tbl$p) & tbl$p < cor_p_max

  # Level 3 / species-level gate: species-wide pooled fit, the final fallback for every
  # individual, and (no minimum sample size beyond the species-inclusion filter already
  # applied upstream) the reliability check deciding whether this species gets
  # SLI-estimated at all.
  pooled <- fit_group(d)
  if (!isTRUE(reliable(pooled, n_min = 0))) {
    return(df %>% dplyr::mutate("{rename_col}" := NA_real_))
  }
  if (length(covariates) == 0) {
    return(df %>% dplyr::mutate("{rename_col}" := {{ Append }} * (M0 / {{ Mass }})^pooled$slope))
  }

  # Collapse unknown-coded/NA values in each covariate to an explicit "Unk" class.
  for (v in covariates) {
    raw <- as.character(d[[v]])
    d[[paste0(".", v, "_cls")]] <- dplyr::if_else(is.na(raw) | raw %in% unknown_codes, "Unk", raw)
  }
  cls_cols <- paste0(".", covariates, "_cls")

  # Level 1: combined-cell fit, only attempted with two covariates.
  if (length(covariates) == 2) {
    cell_tbl <- d %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cls_cols))) %>%
      dplyr::group_modify(~ fit_group(.x)) %>%
      dplyr::ungroup()
    cell_tbl$l1_pass  <- reliable(cell_tbl, n_min_cell)
    cell_tbl$l1_slope <- cell_tbl$slope
    d <- d %>% dplyr::left_join(cell_tbl %>% dplyr::select(dplyr::all_of(cls_cols), l1_slope, l1_pass), by = cls_cols)
  } else {
    d$l1_slope <- NA_real_
    d$l1_pass  <- FALSE
  }

  # Level 2: one marginal fit per covariate, masked to NA wherever it doesn't pass, so a
  # per-row average (below) automatically includes only the passing marginal(s) -- reducing
  # to that one value, unblended, when only one passes.
  for (v in covariates) {
    cls <- paste0(".", v, "_cls")
    marg_tbl <- d %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cls))) %>%
      dplyr::group_modify(~ fit_group(.x)) %>%
      dplyr::ungroup()
    marg_tbl$pass <- reliable(marg_tbl, n_min_marginal)
    marg_tbl[[paste0(".", v, "_marg_masked")]] <- dplyr::if_else(marg_tbl$pass, marg_tbl$slope, NA_real_)
    d <- d %>% dplyr::left_join(marg_tbl %>% dplyr::select(dplyr::all_of(cls), dplyr::ends_with("_marg_masked")), by = cls)
  }

  masked_cols <- paste0(".", covariates, "_marg_masked")
  marg_avg    <- rowMeans(as.data.frame(d[masked_cols]), na.rm = TRUE)  # NaN where none passed
  n_marg_pass <- rowSums(!is.na(as.data.frame(d[masked_cols])))

  b_sli <- dplyr::case_when(
    d$l1_pass %in% TRUE ~ d$l1_slope,
    n_marg_pass >= 1    ~ marg_avg,
    TRUE                ~ pooled$slope
  )

  df %>% dplyr::mutate("{rename_col}" := {{ Append }} * (M0 / {{ Mass }})^b_sli)
}

# calc_lambda function: calculate the empirical coefficients of variation
calc_lambda <- function(x, y){ 
  cv_y <- sd({{ y }}) / mean({{ y }}) 
  cv_x <- sd({{ x }}) / mean({{ x }}) 
  (cv_y^2) / (cv_x^2) 
} 

## Generate correlated mass and wing data
# Helpful for playing around to see how slopes vary with different relationships of mass and wing
# var1 = appendage, mass = mass
### Now a thin wrapper over sliR::sim_correlated(). Raw scale, no log transform and no allometric target, so it is the counterpart to gen_data() for building intuition about how a fitted slope responds to correlation and to error on one trait but not the other.
# sliR names the appendage column Append; rename it back to Appendage so existing call sites are unaffected. n was hardcoded at 3000 and is now an argument, defaulting to that.
gen_cor_vars <- function(r_12, mu_append, mu_mass, sd_append, sd_mass,
                         transient_error_append, transient_error_mass, meas_error,
                         n = 3000) {
  sliR::sim_correlated(
    n         = n,
    r         = r_12,
    mu_append = mu_append,
    mu_mass   = mu_mass,
    sd_append = sd_append,
    sd_mass   = sd_mass,
    meas_error             = meas_error,
    transient_error_append = transient_error_append,
    transient_error_mass   = transient_error_mass
  ) %>%
    dplyr::rename(Appendage = Append)
}

# Classify shapeshifting direction (Bergmann's, Inverse Bergmann's, Mixed, Stable)
# from tidy lm output with mass and wing gradient models per species.
# mass_dir = TRUE means mass decreases along gradient (Bergmann's direction for mass).
# wing_dir = TRUE means wing increases along gradient (Allen's direction for wing).
classify_direction <- function(mods_tbl, p_threshold = 0.05,
                               species_col = "species_",
                               mass_dv = "mass", wing_dv = "wing") {
  mods_tbl %>%
    rename(species_ = !!sym(species_col)) %>%
    mutate(sig = p.value < p_threshold) %>%
    group_by(species_) %>%
    summarise(
      n_sig     = sum(sig),
      mass_dir  = estimate[dv == mass_dv] < 0,
      wing_dir  = estimate[dv == wing_dv] > 0,
      mass_sig  = sig[dv == mass_dv],
      wing_sig  = sig[dv == wing_dv],
      Sig_trait = case_when(
        n_sig == 0 ~ "Neither",
        n_sig == 2 ~ "both",
        TRUE       ~ dv[sig][1]
      ),
      .groups = "drop"
    ) %>%
    mutate(
      sole_decr = case_when(
        mass_sig & !wing_sig ~  mass_dir,
        !mass_sig & wing_sig ~ !wing_dir,
        .default = NA
      ),
      Direction = case_when(
        n_sig == 0                             ~ "Stable",
        n_sig == 1 & sole_decr == TRUE         ~ "Bergmann's",
        n_sig == 1 & sole_decr == FALSE        ~ "Inverse Bergmann's",
        n_sig == 2 &  mass_dir & !wing_dir     ~ "Bergmann's",
        n_sig == 2 & !mass_dir &  wing_dir     ~ "Inverse Bergmann's",
        n_sig == 2 &  mass_dir &  wing_dir     ~ "Mixed - Wingier",
        n_sig == 2 & !mass_dir & !wing_dir     ~ "Mixed - Stouter",
        TRUE ~ "Check"
      )
    ) %>%
    dplyr::select(-sole_decr) %>%
    rename(!!sym(species_col) := species_)
}

# Shared boxplot skeleton for a Parms_tbl4-shaped df (Model, Temp_eff, b_temp_inc, Strength,
# Scaling columns): one panel per Temp_eff category, methods on the x-axis. Used by the main
# text's fig-compare-approaches and by supplementary_info.qmd's ratio-of-logs robustness check
# (same plot, fed the alt-ratio-substituted data) -- moved here once a second call site existed,
# per this project's "extract a shared helper the first time logic is duplicated" convention.
plot_approaches <- function(df, x_txt_size = 9, legend.pos = "top") {
  df <- df %>%
    filter(Temp_eff != "Proportionally larger") %>%
    mutate(Model = str_replace_all(Model, x_labs))

  ggplot(df, aes(x = Model, y = b_temp_inc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = .5, size = 1.9,
                aes(color = Strength, shape = Scaling)) +
    facet_wrap(vars(Temp_eff)) +
    labs(x = NULL, y = expression(hat(beta)[T] ~ "on relative appendage length"),
         color = "Strength", shape = "Scaling") +
    guides(shape = guide_legend(nrow = 1)) +
    theme(
      axis.text.x = element_text(size = x_txt_size, vjust = .58, angle = 55),
      legend.position = legend.pos, legend.box = "vertical"
    )
}
