#' ---
#' title: "Replication: Time Series Growth Curves (`tsgc`)"
#' author: "Ashby, Harvey, Kattuman, Tang, Thamotheram"
#' =======================================================

# ========
# Contents
# ========
# 0. Overview, Model Equations & Data Provenance
#    - One-page results summary
#    - Model equations and assumptions
#    - Data provenance and package pinning
#
# 1. Setup & Utilities
#    - Parameters, libraries, theme, directories, helpers
#
# 2. Baseline Gompertz Model: Gauteng
#    - Data loading, inspection
#    - Gompertz variants (free q, AR(1), fixed q)
#    - Forecasts, holdout accuracy
#
# 3. Gompertz with Exogenous Regressors (xpred)
#    - Weather regressors
#    - Supplying future xpred (idx_series, CSV)
#    - Forecast comparisons
#
# 4. Reproduction Number (R_t)
#    - Mapping Gompertz estimates to R_t
#
# 5. Reinitialisation for Subsequent Waves
#    - Trigger diagnostics (slope & uncertainty)
#    - Forecasts with reinitialisation
#
# 6. Leading Indicator Model: England (Daily)
#    - Baseline leading-indicator estimation
#    - With weather regressors
#
# 7. Comparing Leading Indicator and Gompertz Models (UK-Italy)
#    7.1 Case 1: First peak window (2020-02-25 to 2020-04-01)
#        - UK-only Gompertz vs Italy->UK leading indicator
#    7.2 Cross-validation (same window as Case 1)
#        - Gompertz vs leading-indicator models with lags
#    7.3 Case 2: Extended window (2020-02-25 to 2020-04-15)
#        - Re-comparison under longer sample
#
# 8. Extensions to Other Frequencies
#    8.1 Quarterly: Wii sales (Gompertz)
#    8.2 Quarterly: Wii->Switch (leading indicator)
#    8.3 Monthly: Plus500 (Gompertz)
#    8.4 Monthly: DEGIRO->AvaTrade (leading indicator)
#    8.5 Annual: 3DS (Gompertz)
#    8.6 Annual: Wii->3DS (leading indicator)
#
# 9. Appendix: Controlling the Plot X-Axis with idx_axis_opts()
#    - mode = "date" / "position" / "steps" / "time_since"
#    - info_box, pattern_n
#
# 10. Appendix: Compound and Heterogeneous Calendar Steps
#     10.1 idx_step() / idx_calendar_step(): a single compound step
#     10.2 multi_step_pattern() / idx_calendar_multi_step(): heterogeneous cycles
#     10.3 idx_offset_to_pos(): non-calendar idx_calendar anchors
#
# 9. Limitations, Diagnostics & Exported-File Notes
#    - Illustrative vs retrospective vs operational forecasts
#    - Cross-validation scope and interval calibration
#    - Convergence/residual diagnostics checklist
#    - Exported CSV file labelling
#
# ==========================================
#
# This script works with `idx_series` (integer-position
# data) rather than directly with `xts` objects: estimation and
# forecasting are defined purely in terms of position, independent of
# calendar frequency or gaps. Calendar time is reintroduced only for
# plotting, via an optional `idx_calendar` that is attached to each
# model with the `calendar` argument. `xts_to_idx()` converts a
# calendar-indexed `xts`/`zoo` object (e.g. the package's bundled
# datasets) to this representation in one step, returning both the
# `idx_series` and a matching `idx_calendar`; `idx_to_pos()` translates
# a calendar date into the integer position that `start`/`end`/
# `reinit.idx`/`n.lag` expect.

#' 
#' # 0. Overview, Model Equations & Data Provenance
#' 
#' ## 0.1 One-page results summary
#' 
#' | Example | Frequency | Model | Estimation window | Horizon |
#' |---|---|---|---|---|
#' | Gauteng cases (free q / fixed q=0.005) | Daily | Dynamic Gompertz | 2021-02-01 to 2021-05-03 | 14d |
#' | Gauteng cases + weather (q=0.005, fixed to match holdout) | Daily | Dynamic Gompertz + xpred | 2021-02-01 to 2021-05-03 | 14d |
#' | England hospital admissions (+ weather) | Daily | Leading indicator (+ xpred) | 2021-04-30 to 2021-07-24 | 14d |
#' | UK-Italy Case 1 / Case 2 | Daily | Gompertz vs leading indicator | 2020-02-25 to 2020-04-01 / 04-15 | 14d |
#' | Wii / Plus500 / DEGIRO-AvaTrade | Quarterly / Monthly | Gompertz / leading indicator | see relevant section | 4 |
#' | 3DS / Wii->3DS (annual, Q4-dated) | Annual | Gompertz / leading indicator | 2011-12 to 2018-12 | 2y |
#' 
#' Fill in forecast/accuracy figures for each row from this run's own
#' output (`results/Tables`); figures are intentionally omitted here since
#' they depend on run-specific data vintage. Read this table alongside the
#' limitations in Section 9 -- several rows (the annual examples in
#' particular) are illustrative, not empirically validated, forecasts.
#' 
#' ## 0.2 Model equations and assumptions
#' 
#' **Dynamic Gompertz model** (Sections 2-3, 5, 8.1, 8.3, 8.5): models the
#' log-growth rate of a cumulative series as a locally-linear trend in a
#' state-space form, following Harvey & Kattuman (2021), "A farewell to R:
#' time-series models for tracking and forecasting epidemics", Journal of
#' the Royal Society Interface, 18. <http://doi.org/10.1098/rsif.2021.0179>.
#' 
#' - `q` is the signal-to-noise ratio governing how much the slope state
#'   drifts period to period. It can be estimated freely or fixed; **the
#'   value used must be reported alongside every table or forecast**,
#'   since free-q and fixed-q fits are different models (see the corrected
#'   `model_weather` specification in Section 3, which now fixes
#'   `q = q.default` to match its paired holdout model).
#' - `xpred` regressors enter as an `SSMregression` component. When
#'   supplying future values via `supply_xpred.new()`, the model performs
#'   a **date match** between the regressor series' index and the forecast
#'   origin (`end.date`), not a positional/row-number match. Confirm index
#'   classes agree (e.g. `Date` vs `POSIXct`) between the model's own
#'   dates and the supplied regressor series.
#' 
#' **Leading-indicator model** (Sections 6-7, 8.2, 8.4, 8.6): a bivariate
#' extension in which a lead series' lagged values inform the target
#' series' state, with the lag length `n.lag` specified by the user (see
#' Section 7.2 on lag selection, and Section 9 on its validation scope).
#' 
#' **Reproduction number, R_t** (Section 4): computed as
#' \(R_t = \exp(g_{y,t} \cdot \text{gen\_int})\), where \(g_{y,t}\) is the
#' fitted daily log-growth rate and `gen_int` (default 4 days) is an
#' **assumed**, not estimated, generation interval. `ndays` controls how
#' many of the most recent daily R_t values are returned/plotted -- it is
#' a truncation parameter, not a smoothing window. See the `gen_int`
#' sensitivity table added in Section 4.2.
#' 
#' ## 0.3 Data provenance and package pinning
#' 
#' All series used below (`gauteng`, `gauteng_weather_2021`, `england`,
#' `england_weather_2021`, `nintendo_sales`, `etrading_apps`, and the
#' Italy comparison series) are built-in datasets shipped with the `tsgc`
#' package itself, used as-is except for the column subsetting and
#' frequency conversions documented at each point of use (e.g. the
#' Q4/year-end conversion of `nintendo_sales` in Section 8.5).
#' 
#' This document requires `tsgc 2.0.0`. For a fully pinned replication,
#' record here: the exact repository URL and commit SHA used, the
#' installation command (e.g.
#' `remotes::install_github("<org>/tsgc", ref = "<commit-sha>")`), and an
#' `renv.lock` capturing the full dependency graph (`KFAS`, `dplyr`,
#' `xts`, etc.), not just `tsgc` itself.
#' 

#' 
#' # 1. Setup & Utilities
#' 
#' Centralise parameters, load libraries, set global options, define paths,
#' and provide helper functions for saving plots and controlling chunk defaults.
#' 

## ---- 1-0-preamble, include=TRUE---
# Start the reproducible analysis script. This preamble establishes shared settings before any model is estimated.

# ====================
# 1. Setup & Utilities
# ====================

# ---- 1.1 Parameters & Toggles ----
# Define global switches and defaults so forecast horizons, plot sizes, confidence levels, and output behaviour can be changed in one place.
# These booleans control whether figures and tables are written to disk; set them to FALSE for a dry run.
# NOTE: validate_saved_figure() (see save_plot()/safe_ggsave() below) only
# ever runs when a figure file is actually written, i.e. only when
# SAVE_PLOTS is TRUE. With the default FALSE below, no PNGs are written
# and the render-validation check is correspondingly not exercised. Set
# SAVE_PLOTS <- TRUE for any run whose purpose is to confirm figures
# actually render (e.g. before publishing or archiving this document).
SAVE_PLOTS   <- FALSE
SAVE_TABLES  <- FALSE
FIG_WIDTH    <- 10
FIG_HEIGHT   <- 7
FIG_DPI      <- 300
CONF_LEVEL   <- 0.68  # Coverage proportion for plotted/exported uncertainty intervals.
# 0.68: ~one-standard-error interval under Gaussian assumptions.

# Core analysis parameters 
# The default 14-step horizon is used repeatedly for the daily examples unless a later section overrides it.
n.forecasts.default <- 14
q.default           <- 0.005
plt.length.default  <- 30

# Reproduction Number (R_t) parameters
# These epidemiological assumptions are used only when translating growth dynamics into R_t.
gen_int <- 4   # assumed generation interval, how many days between infection in one case 
# and secondary infections.
ndays   <- 7   # Number of days used to smooth/aggregate the growth signal when estimating R_t.
# A 7-day window helps reduce daily reporting noise in COVID case data.


# Default estimation window used in Section 2, as calendar dates. These
# are translated to integer positions via idx_to_pos() once the Gauteng
# calendar (gauteng_cal) has been constructed in Section 2.1.
est.start.1.date <- as.Date("2021-02-01")
est.end.1.date   <- as.Date("2021-05-03")

# User note: Sections 1.2 to 1.7 contain setup code for loading packages,
# defining output folders, saving plots/tables, exporting CSV files, and setting
# knitr options. Readers interested mainly in the modelling examples can skip
# ahead to Section 2 after these setup commands have been run.

# ---- 1.2 Libraries (quiet require) ----
# Load all package dependencies quietly. The helper keeps the console output focused while attaching the packages needed below.
# Define a small loader so every required package is attached with suppressed startup messages.
safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install before running.", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# List all packages used by the script;  missing dependencies are easier to diagnose.
libs <- c(
  "tsgc","KFAS","dplyr","ggplot2","ggthemes",
  "zoo","xts","gridExtra","here",
  "tidyr","abind","scales","grid","png"
)
invisible(lapply(libs, safe_library))

# ---- 1.2b Reproducibility: pin and verify the installed tsgc build ----
# Naming a commit SHA in prose is not the same as verifying it against
# what is actually installed, and a full test-suite/R-CMD-check run
# does not by itself pin dependency versions for anyone re-running this
# script later. This block: (a) records the exact installed tsgc
# version/commit actually used for this run, so results can be tied to
# a specific build rather than "whatever tsgc happened to be installed
# on the day"; (b) writes a full sessionInfo() to disk alongside any
# other outputs; and (c) attempts an renv.lock snapshot if the renv
# package is available, so dependency versions (not just tsgc) are
# pinned too. None of this was executed against the actual installed
# package in this environment (no R available here) - this only wires
# up the pinning/verification infrastructure; the printed
# version/commit and the written renv.lock must be checked/generated by
# actually running this block.
EXPECTED_TSGC_COMMIT <- NULL  # set to the exact commit SHA this replication was built against, e.g. "a1b2c3d4..."

tsgc_desc <- tryCatch(utils::packageDescription("tsgc"), error = function(e) NULL)
tsgc_version <- if (!is.null(tsgc_desc)) tsgc_desc$Version else NA_character_
`%||%` <- function(x, y) if (is.null(x)) y else x
# GithubSHA1/RemoteSha are populated by devtools::install_github()/
# remotes::install_github(); a CRAN or local install won't have these,
# in which case only the package Version is available to check.
tsgc_commit <- if (!is.null(tsgc_desc)) {
  tsgc_desc$GithubSHA1 %||% tsgc_desc$RemoteSha %||% NA_character_
} else NA_character_

message("Installed tsgc version: ", tsgc_version)
message("Installed tsgc commit/remote SHA: ", if (is.na(tsgc_commit)) "unavailable (not installed via install_github()/install_git())" else tsgc_commit)

if (!is.null(EXPECTED_TSGC_COMMIT)) {
  if (is.na(tsgc_commit) || !identical(tsgc_commit, EXPECTED_TSGC_COMMIT)) {
    warning(
      "Installed tsgc commit ('", tsgc_commit, "') does not match the ",
      "commit this replication was pinned to ('", EXPECTED_TSGC_COMMIT,
      "'). Results below may not match the original replication exactly."
    )
  } else {
    message("Installed tsgc commit matches EXPECTED_TSGC_COMMIT.")
  }
}

# NOTE: the sessionInfo()/renv.lock pinning step itself (using
# ensure_dir(), results_dir, and base_path) runs later, in Section 1.4,
# once those are actually defined - see below.

# ---- 1.3 Global Options & Theme ----
# Set plotting and printing defaults so figures are visually consistent across the replication.
theme_set(ggthemes::theme_economist_white(gray_bg = FALSE, base_size = 16))
options(scipen = 7)

# ---- 1.4 Paths & Directories ----
# Construct local output folders relative to the current working directory. Avoids hard-coded user-specific paths.
# Set the working directory to the package root with 
# Session > Set Working Directory > Choose Directory...
# This ensures results/ is created in the project root 
base_path   <- getwd()
results_dir <- file.path(base_path, "results")
tables_dir  <- file.path(results_dir, "Tables")
images_dir  <- file.path(results_dir, "Images")

# Create output directories if they do not already exist.
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- tryCatch({
      dir.create(path, recursive = TRUE); TRUE
    },
    error = function(e) {
      message("Failed to create: ", path, " :: ", e$message); FALSE
    })
    if (!ok) stop("Could not create required directory: ", path)
  }
}
invisible(lapply(list(results_dir, tables_dir, images_dir), ensure_dir))

# ---- 1.4b Reproducibility: session/dependency snapshot ----
# Runs after ensure_dir(), results_dir, and base_path are defined.
if (SAVE_TABLES) {
  ensure_dir(results_dir)
  writeLines(
    capture.output(sessionInfo()),
    con = file.path(results_dir, "sessionInfo.txt")
  )
  message("Saved sessionInfo.txt (full package/version manifest for this run).")
  
  if (requireNamespace("renv", quietly = TRUE)) {
    tryCatch({
      renv::snapshot(project = base_path, prompt = FALSE)
      message("Saved renv.lock via renv::snapshot().")
    }, error = function(e) {
      message("renv::snapshot() failed: ", conditionMessage(e),
              " - install/configure renv to pin dependency versions.")
    })
  } else {
    message("renv not installed: skipping renv.lock snapshot. Install ",
            "renv and re-run with SAVE_TABLES = TRUE to pin dependency versions.")
  }
}

# ---- 1.5 Plot/Save Helpers ----
# Define wrappers for saving plots. These centralise image dimensions, resolution, and the SAVE_PLOTS toggle.
# Save a ggplot with common width, height, and DPI settings when SAVE_PLOTS is TRUE.
safe_ggsave <- function(plot, filename, width = FIG_WIDTH, 
                        height = FIG_HEIGHT, dpi = FIG_DPI) {
  # Returns FALSE (no file written) when SAVE_PLOTS is off, and TRUE
  # once ggsave() has completed - this return value, not a re-check of
  # SAVE_PLOTS, is what the caller uses to decide whether to validate.
  if (!SAVE_PLOTS) return(invisible(FALSE))
  ggplot2::ggsave(filename = filename, plot = plot, 
                  width = width, height = height, dpi = dpi)
  message("Saved plot: ", normalizePath(filename, winslash = "/"))
  invisible(TRUE)
}

# Print each plot in interactive/rendered output and optionally save it to the Images folder.
# Validation is tied to whether a file was actually written (i.e. to
# SAVE_PLOTS via safe_ggsave's own no-op), not to a second, separately
# maintained SAVE_PLOTS check here. With the default SAVE_PLOTS = FALSE,
# no file is written and validate_saved_figure() is correctly skipped;
# whenever SAVE_PLOTS = TRUE, every saved figure is validated
# unconditionally, so a zero-byte figure is always caught, not only
# when a second flag happens to also be TRUE.
save_plot <- function(p, fname = NULL) {
  print(p)
  if (!is.null(fname) && inherits(p, "ggplot")) {
    filepath <- file.path(images_dir, fname)
    saved <- safe_ggsave(p, filepath)
    if (isTRUE(saved)) validate_saved_figure(filepath)
  }
  invisible(p)
}

# Validate that a saved figure actually rendered: the file must exist,
# be non-empty, decode as an image, and exceed a plausible minimum
# pixel size. A positive file size alone is not sufficient to confirm
# a figure rendered correctly; the file must also decode.
validate_saved_figure <- function(filepath, min_dim = 50) {
  if (!file.exists(filepath)) {
    stop("Figure was not created: ", filepath)
  }
  sz <- file.info(filepath)$size
  if (is.na(sz) || sz <= 0) {
    stop("Figure is zero bytes (failed render): ", filepath)
  }
  dims <- tryCatch(
    {
      raster <- png::readPNG(filepath, info = TRUE)
      dim(raster)[1:2]
    },
    error = function(e) {
      stop("Figure does not decode as a valid PNG: ", filepath,
           " (", conditionMessage(e), ")")
    }
  )
  if (any(is.na(dims)) || any(dims < min_dim)) {
    stop("Figure decodes but is implausibly small (", 
         paste(dims, collapse = "x"), " px): ", filepath)
  }
  invisible(TRUE)
}

# ---- 1.5b Model Diagnostics Helper ----
# print_model_diagnostics() lives in the package (utils.R). It operates
# purely on public FilterResults/FilterResultsLI accessors and is
# generically useful, not specific to this replication exercise. It
# reports the log-likelihood and recursive-residual
# diagnostics - none of which summary() on these classes reports.

# ---- 1.6 CSV Export Helpers ----
# Define CSV-export helpers for forecasts, filtered states, growth rates, R_t outputs, and the manifest used to document exported files.
# Since estimation works on idx_series (integer positions) rather
# than xts, and models carry an optional idx_calendar, dates for export
# are recovered via idx_to_date(calendar, idx_positions(x)) when a
# calendar is available, falling back to raw integer positions otherwise.
format_csv_dates <- function(calendar, positions) {
  if (!is.null(calendar)) {
    format(idx_to_date(calendar, positions), "%Y-%m-%d")
  } else {
    as.character(positions)
  }
}

# Turn the confidence level into a suffix, for example 68pct, for unambiguous exported column names.
confidence_suffix <- function(confidence.level) {
  as.character(round(confidence.level * 100))
}

# Export an idx_series object with its positions restored as an explicit Date (or position) column.
write_idx_csv <- function(x, calendar, file, columns) {
  out <- data.frame(
    Date = format_csv_dates(calendar, idx_positions(x)),
    as.matrix(idx_values(x)),
    check.names = FALSE
  )
  names(out) <- c("Date", columns)
  # Rows with an exact zero standard error (typically diffuse
  # initialisation at the start of the sample) read as false certainty
  # if exported as-is - the associated confidence bounds collapse to
  # the point estimate rather than reflecting genuine (unbounded/
  # undefined) diffuse-state uncertainty. Rather than only flagging
  # these rows, null out the SE and any interval bounds derived from
  # it, so a false-precision zero can't be silently consumed
  # downstream (e.g. averaged into a coverage or width statistic).
  # diffuse_flag is still kept, recording which rows were affected.
  err_col <- grep("std_error", names(out), value = TRUE)
  if (length(err_col) == 1 && any(out[[err_col]] == 0, na.rm = TRUE)) {
    diffuse_rows <- out[[err_col]] == 0
    out$diffuse_flag <- diffuse_rows
    bound_cols <- grep("_lower_|_upper_|^lower$|^upper$", names(out), value = TRUE)
    out[diffuse_rows, err_col] <- NA_real_
    if (length(bound_cols)) {
      out[diffuse_rows, bound_cols] <- NA_real_
    }
    message("Note: ", sum(diffuse_rows), " row(s) in ", basename(file),
            " had zero standard error (diffuse initialisation); SE and any",
            " derived interval bounds were set to NA and diffuse_flag was set.")
  }
  write.csv(out, file = file, row.names = FALSE)
  message("Saved table: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(out)
}

# Export the main forecast, filtered states, growth-rate estimates, and uncertainty intervals for a fitted model.
write_results_clear <- function(res, res.dir, n.ahead, model_slug, target_slug,
                                confidence.level = CONF_LEVEL) {
  ensure_dir(res.dir)
  
  ci <- confidence_suffix(confidence.level)
  forecast_col <- paste0("forecast_", target_slug)
  forecast_lower_col <- paste0("forecast_lower_", ci)
  forecast_upper_col <- paste0("forecast_upper_", ci)
  growth_lower_col <- paste0("growth_rate_lower_", ci)
  growth_upper_col <- paste0("growth_rate_upper_", ci)
  
  calendar <- res$calendar
  
  y.hat.diff <- res$predict_level(
    n.ahead = n.ahead,
    confidence.level = confidence.level,
    sea.on = TRUE
  )
  write_idx_csv(
    y.hat.diff,
    calendar,
    file.path(res.dir, paste0(model_slug, "_", target_slug, "_forecast.csv")),
    c(forecast_col, forecast_lower_col, forecast_upper_col)
  )
  
  y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
  filtered.level <- y.hat.all$level.t.t
  filtered.slope <- y.hat.all$slope.t.t
  a.t.t <- y.hat.all$a.t.t
  P.t.t <- y.hat.all$P.t.t
  
  idx.level <- grep("level", colnames(a.t.t))[1]
  idx.slope <- grep("slope", colnames(a.t.t))[1]
  delta.var  <- P.t.t[idx.level, idx.level, ]
  gamma.var  <- P.t.t[idx.slope, idx.slope, ]
  delta.gamma.cov <- P.t.t[idx.level, idx.slope, ]
  delta.std.err <- sqrt(delta.var)
  gamma.std.err <- sqrt(gamma.var)
  
  # The *_filtered.csv files should contain only in-sample filtered
  # states, not rows extending into the forecast period. Restrict to
  # positions at or before the model's own fitted sample end.
  filtered.end.pos <- tail(res$index, 1)
  
  delta <- idx_series(
    cbind(
      delta_log_growth_level = idx_values(filtered.level),
      delta_std_error = as.numeric(delta.std.err)
    ),
    start = filtered.level$start
  )
  write_idx_csv(
    get_timeframe(delta, filtered.level$start, filtered.end.pos),
    calendar,
    file.path(res.dir, paste0(model_slug, "_delta_filtered.csv")),
    c("delta_log_growth_level", "delta_std_error")
  )
  
  gamma <- idx_series(
    cbind(
      gamma_trend_slope = idx_values(filtered.slope),
      gamma_std_error = as.numeric(gamma.std.err)
    ),
    start = filtered.slope$start
  )
  write_idx_csv(
    get_timeframe(gamma, filtered.slope$start, filtered.end.pos),
    calendar,
    file.path(res.dir, paste0(model_slug, "_gamma_filtered.csv")),
    c("gamma_trend_slope", "gamma_std_error")
  )
  
  # g_t = exp(delta_t) + gamma_t. Following Harvey & Kattuman (2021,
  # Sec 3.1), the contribution of Var(delta_t) to Var(g_t) is negligible
  # relative to Var(gamma_t) once the epidemic is underway, so the
  # sampling variability of g_t is taken as that of gamma_t alone:
  #   Var(g_t) ~= Var(gamma_t)
  e.delta <- exp(idx_values(filtered.level))
  fitted.growth <- e.delta + idx_values(filtered.slope)
  growth.var <- as.numeric(gamma.var)
  growth.se <- sqrt(pmax(growth.var, 0))
  ci.offset <- stats::qnorm((1 - confidence.level) / 2) *
    growth.se %o% c(1, -1)
  growth.ci <- idx_series(
    cbind(
      fitted_incidence_growth_rate = as.numeric(fitted.growth),
      lower = as.numeric(as.numeric(fitted.growth) + ci.offset[, 1]),
      upper = as.numeric(as.numeric(fitted.growth) + ci.offset[, 2])
    ),
    start = filtered.level$start
  )
  write_idx_csv(
    growth.ci,
    calendar,
    file.path(res.dir, paste0(model_slug, "_", target_slug, "_growth_rate.csv")),
    c("fitted_incidence_growth_rate", growth_lower_col, growth_upper_col)
  )
  
  invisible(TRUE)
}

# Write a manifest describing each CSV output so downstream users know what each file contains.
write_csv_manifest <- function(res.dir) {
  ci <- confidence_suffix(CONF_LEVEL)
  manifest <- data.frame(
    file = c(
      "gauteng_gompertz_q005_new_cases_forecast.csv",
      "gauteng_gompertz_q005_new_cases_growth_rate.csv",
      "gauteng_gompertz_q005_delta_filtered.csv",
      "gauteng_gompertz_q005_gamma_filtered.csv",
      "gauteng_gompertz_q005_rt.csv",
      "england_leading_indicator_hospital_admissions_forecast.csv",
      "england_leading_indicator_hospital_admissions_growth_rate.csv",
      "england_leading_indicator_delta_filtered.csv",
      "england_leading_indicator_gamma_filtered.csv"
    ),
    description = c(
      "Gauteng 14-day forecast of daily new cases, with prediction interval bounds.",
      "Gauteng fitted incidence growth rate, with confidence interval bounds.",
      "Gauteng filtered delta state: log-growth level and standard error.",
      "Gauteng filtered gamma state: trend slope and standard error.",
      "Gauteng effective reproduction number estimate, with confidence interval bounds.",
      "England 14-day forecast of hospital admissions, with prediction interval bounds.",
      "England fitted hospital-admission growth rate, with confidence interval bounds.",
      "England filtered delta state: log-growth level and standard error.",
      "England filtered gamma state: trend slope and standard error."
    ),
    columns = c(
      paste0("Date, forecast_new_cases, forecast_lower_", ci,
             ", forecast_upper_", ci),
      paste0("Date, fitted_incidence_growth_rate, growth_rate_lower_", ci,
             ", growth_rate_upper_", ci),
      "Date, delta_log_growth_level, delta_std_error",
      "Date, gamma_trend_slope, gamma_std_error",
      paste0("Date, Rt, Rt_lower_", ci, ", Rt_upper_", ci),
      paste0("Date, forecast_hospital_admissions, forecast_lower_", ci,
             ", forecast_upper_", ci),
      paste0("Date, fitted_incidence_growth_rate, growth_rate_lower_", ci,
             ", growth_rate_upper_", ci),
      "Date, delta_log_growth_level, delta_std_error",
      "Date, gamma_trend_slope, gamma_std_error"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, file.path(res.dir, "csv_manifest.csv"), row.names = FALSE)
  message("Saved table: ",
          normalizePath(file.path(res.dir, "csv_manifest.csv"),
                        winslash = "/", mustWork = FALSE))
  invisible(manifest)
}

# Convenience for plot windows expressed as integer positions
# Helper for plot windows: move back k positions from the last available position.
tail_pos_minus <- function(pos_vec, k)
  if (length(pos_vec)) tail(pos_vec, 1) - k else NA

# Takes the last position and goes back k positions;
# if the vector is empty, returns NA.

# ---- 1.7 knitr Defaults ----
# Configure knitr chunk defaults when the script is rendered as a vignette or report.
# NOTE: warnings are shown (not globally suppressed). The prior
# `warning = FALSE` hid many "rows removed" notices (some 40-50 rows)
# that should be understood and handled locally rather than silenced.
knitr::opts_chunk$set(
  echo       = TRUE,
  message    = TRUE,
  warning    = TRUE,
  fig.align  = "center",
  fig.width  = FIG_WIDTH,
  fig.height = FIG_HEIGHT
)

# ===================================
# 2. Baseline Gompertz Model: Gauteng
# ===================================

#' # 2. Baseline Gompertz Model: Gauteng
#' 
#' The data used in this replication relate to daily COVID-19 cases in Gauteng
#' province in South Africa. We begin by fitting a baseline dynamic Gompertz model
#' to cumulative COVID-19 cases. We then estimate variants of the baseline model,
#' generate forecasts, and evaluate forecast accuracy.
#' 
#' ## 2.1 Data
#' 
#' Load the Gauteng dataset and convert it to idx_series/idx_calendar form
#' for modelling.
#' 

## ---- 2.1 Data -----------------
# Load the built-in Gauteng dataset from tsgc (ships as an xts object).
data(gauteng, package = "tsgc")

# Convert to idx_series (integer-position data) plus a matching
# idx_calendar that records how positions map back to real dates. All
# estimation/forecasting below works on gauteng_idx; gauteng_cal is
# carried through models purely so that plots can be labelled with dates.
conv <- xts_to_idx(gauteng)
gauteng_idx <- conv$series
gauteng_cal <- conv$calendar
gauteng_cal$anchor_name <- "first recorded case"

# Cumulative cases is the response series for the baseline Gompertz model.
cumulative_cases <- gauteng_idx

# Translate the calendar-date estimation window (Section 1.1) into
# integer positions now that gauteng_cal is available.
est.start.1 <- idx_to_pos(gauteng_cal, est.start.1.date)
est.end.1   <- idx_to_pos(gauteng_cal, est.end.1.date)

#' 
#' ## 2.2 Quick Inspection of Data
#' 
#' Visual sanity check of the series and index.
#' 

## ---- 2.2 Data Inspection ------------------------------------------
# Build a quick model object and plot the raw series as a visual check before formal estimation.
# Construct a model object for plotting the series before formal estimation.
mod1 <- tsgc::SSModelDynamicGompertz(Y = cumulative_cases, calendar = gauteng_cal)
# Plot the input data to check scale, timing, and obvious data problems.
p <- plot(mod1, title = "Gauteng daily cases", series.name = "Cases")
print(p)
save_plot(p, "gauteng_cases_MA.png")

#' 
#' ## 2.3 Estimation: Gompertz Model Variants
#' 
#' Three Gompertz model variants:
#' 
#' - (a) Diffuse / free q: flexible baseline
#' - (b) AR(1) slope: smoother slope dynamics
#' - (c) Fixed q: chosen based on experience
#' 

## ---- 2.3 Estimation ---------------------------
# Estimate three dynamic Gompertz variants to compare flexible, AR(1), and fixed-q slope dynamics.
# 2.3a Diffuse prior (free q)
# Specify a dynamic Gompertz model with q estimated from the data, giving a flexible baseline.
model_free <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, start = est.start.1, 
  end = est.end.1, calendar = gauteng_cal
)
# Estimate the free-q model and print its summary for parameter/state diagnostics.
res_free <- tsgc::estimate(model_free); summary(res_free); tsgc::print_model_diagnostics(res_free)

# 2.3b Diffuse prior with AR(1)
# Specify an AR(1) slope variant to impose smoother slope evolution.
model_ar1 <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, ar1 = TRUE, start = est.start.1, 
  end = est.end.1, calendar = gauteng_cal
)
# Estimate the AR(1) variant and inspect the summary.
res_ar1 <- tsgc::estimate(model_ar1); summary(res_ar1); tsgc::print_model_diagnostics(res_ar1)

# 2.3c Fixed q
# Specify the preferred fixed-q model using the shared q.default value.
model_q <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, q = q.default, start = est.start.1, 
  end = est.end.1, calendar = gauteng_cal
)
# Estimate the fixed-q model used for the baseline forecasts.
res_q <- tsgc::estimate(model_q); summary(res_q); tsgc::print_model_diagnostics(res_q)

#' 
#' ## 2.4 Forecasts & Accuracy
#' 
#' Produce forecasts (log growth and levels) from the fixed-q model
#' and evaluate holdout accuracy.
#' 

## ---- 2.4 Forecasts & Accuracy ----
# Forecast from the selected fixed-q model and evaluate a two-week holdout period.
# Reset the forecast horizon and plotting window to the daily defaults for this section.
n.forecasts <- n.forecasts.default
plt.length  <- plt.length.default

# 2.4a Log growth forecast (fixed q)
# Forecast the latent/log growth component from the fixed-q model.
p <- tsgc::plot_log_forecast(
  res_q, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start = tail_pos_minus(res_q$index, plt.length),
  title = "Forecast of log growth rate of cases\n14-days (Gauteng)"
); print(p)

# 2.4b New cases forecast
# Forecast daily new cases in the original data scale with prediction intervals.
p <- tsgc::plot_forecast(
  res_q, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start = tail_pos_minus(res_q$index, plt.length),
  title = "Forecast of new cases\n14-days (Gauteng)", 
  series.name = "Cases"
); print(p)

# 2.4c Holdout accuracy: two weeks prior to end of sample
# Define a holdout estimation end 14 positions before est.end.1
# Define a truncated estimation end position so the final 14 positions can be held out for validation.
est.end.holdout <- est.end.1 - n.forecasts

# Refit ONLY for holdout evaluation on the truncated window
# Refit the fixed-q Gompertz model using only data available before the holdout period.
model_q_holdout <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases,  
  q = q.default, 
  start = est.start.1,
  end   = est.end.holdout,
  calendar = gauteng_cal
)
# Estimate the holdout model before comparing its forecasts with the withheld observations.
res_q_holdout <- tsgc::estimate(model_q_holdout); 
summary(res_q_holdout)
tsgc::print_model_diagnostics(res_q_holdout)

# 2.4c Holdout accuracy plot
# Plot holdout accuracy by comparing forecasts with observed values after the truncated end position.
p <- tsgc::plot_holdout(
  res_q_holdout, Y = cumulative_cases, n.ahead = n.forecasts, 
  confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases\n14-days (Gauteng)", 
  series.name = "Cases"
); print(p)

if (SAVE_TABLES) {
  # Export the baseline Gauteng forecast and filtered-state results as documented CSV files.
  write_results_clear(
    res = res_q,
    res.dir = tables_dir,
    n.ahead = n.forecasts,
    model_slug = "gauteng_gompertz_q005",
    target_slug = "new_cases",
    confidence.level = CONF_LEVEL
  )
  message("Saved clear CSV results for: gauteng_gompertz_q005")
}

#' 
#' Note: The symmetric Mean Absolute Percentage Error (sMAPE) is a scale-free,
#' symmetric accuracy measure that ranges from 0% to 100%. It complements MAPE,
#' which tends to overstate forecast errors when actual values are small.
#' However, sMAPE can also become unstable when both actual and forecast values
#' are very small.
#' 

# ====================================================
# 3. Gauteng: Gompertz Model with Exogenous Regressors
# ====================================================
# 
# Augment the Gompertz model with weather regressors; re-estimate and compare
# forecasts and accuracy.
#

## ---- 3.1 Estimation with Regressors: xpred ----

# In this section we show how to use exogenous variables (xpred),
# such as future weather values, as regressors in the tsgc model.
# This allows us to generate out-of-sample forecasts conditional
# on a specified path for the regressors.

# ------------------------------------------------------------
# 3.1.1 Use built-in weather data as regressors for estimation
# ------------------------------------------------------------

# The example data `gauteng_weather_2021` (included in the tsgc package)
# contain daily weather variables for 2021. We will use:
#   - column 1: Wind speed
#   - column 3: Mean daily temperature
# as regressors over the *same* estimation window as before.

# Load built-in daily Gauteng weather data and convert to idx_series,
# anchored on the Gauteng calendar so that positions line up with
# cumulative_cases.
data(gauteng_weather_2021, package = "tsgc")
gauteng_weather_idx <- xts_to_idx(
  gauteng_weather_2021[, c(1, 3)],
  start.pos = idx_to_pos(gauteng_cal, zoo::index(gauteng_weather_2021)[1])
)$series

# Subset to the estimation window [est.start.1, est.end.1]
# Restrict weather regressors to the same estimation window as the response series.
gauteng_weather_est <- get_timeframe(
  gauteng_weather_idx,
  est.start.1,
  est.end.1
)
head(gauteng_weather_est)

# Fit a Dynamic Gompertz model with weather regressors
# Fit the Gompertz model with weather regressors included through xpred.
# NB. q is fixed at q.default here so that this full-sample specification
# matches the holdout specification (model_q_xpred_holdout, Section 3.1.4),
# which also fixes q = q.default. Estimating q freely in one and fixing it
# in the other would make the forecast and the holdout-accuracy table
# describe two different models rather than the same model evaluated two
# ways. If a free-q comparison is wanted, add it as a clearly separate,
# additionally-labelled specification rather than substituting it here.
model_weather <- tsgc::SSModelDynamicGompertz(
  Y          = cumulative_cases,
  xpred      = gauteng_weather_est,
  q          = q.default,
  start      = est.start.1,
  end        = est.end.1,
  calendar   = gauteng_cal
)

# Estimate the weather-augmented model and inspect its fitted output.
res_weather <- tsgc::estimate(model_weather)
summary(res_weather)
tsgc::print_model_diagnostics(res_weather)

# ---------------------------------------------------
# 3.1.2 Supplying future xpred values for forecasting
# ---------------------------------------------------

# The object `res_weather` contains parameter estimates based on
# *in-sample* weather data (the estimation period). To produce
# out-of-sample forecasts, the model needs xpred values BEYOND
# the estimation window.

# These future xpred values can be:
#   - actual weather forecasts from a provider (e.g. Met Office), or
#   - user-specified scenarios (best/worst-case paths).

# *** ORACLE-FORECAST CAVEAT ***
# The example below instead takes REALISED (observed after the fact)
# weather from `gauteng_weather_idx` for the forecast horizon, not an
# archived weather forecast that would actually have been available at
# the forecast origin. This makes the resulting accuracy figures a
# conditional ("oracle") forecast exercise, not an operational
# validation: it shows what the model would have predicted GIVEN
# perfect future weather knowledge, not what it would have predicted
# using information available at the time. Any accuracy claims below
# should be read with this caveat; a genuine operational evaluation
# would require archived point/ensemble forecasts as of each origin
# date, not realised values.

# Example: take future rows from `gauteng_weather_idx` for the
# forecast horizon of length n.forecasts:
# Extract the future weather path needed for out-of-sample forecasts.
gauteng_weather_future <- get_timeframe(
  gauteng_weather_idx,
  est.end.1 + 1,
  est.end.1 + n.forecasts
)

# Supply these future regressors to the fitted model. The 
# xpred.new field is set directly on
# the fitted FilterResults object.
res_weather$xpred.new <- gauteng_weather_future

# ----------------------------------------------------------
# 3.1.3 Example: reading future xpred values from a CSV file
# ----------------------------------------------------------

# In practice, future xpred values will often be read from a CSV file.
# The CSV must contain:
#   - a column named 'Date'
#   - one or more numeric columns with length n.forecasts.
#
# Here we illustrate with a small inline CSV, 
# read into `gauteng_weather_future_csv`.

# This inline CSV block demonstrates an alternative way to supply future xpred values from a file-like source.
txt <- "
Date,windspd_mtrs_p_sec,temperature_C
2021-05-04,2.03,14.61
2021-05-05,1.54,15.51
2021-05-06,1.94,16.42
2021-05-07,2.38,15.54
2021-05-08,2.57,14.18
2021-05-09,2.65,13.55
2021-05-10,2.19,14.84
2021-05-11,2.08,15.55
2021-05-12,2.07,15.97
2021-05-13,1.92,15.86
2021-05-14,1.85,15.79
2021-05-15,2.86,15.57
2021-05-16,3.46,16.42
2021-05-17,2.39,13.53
"

# Read the demonstration CSV text into a data frame with dates preserved.
gauteng_weather_future_csv <- read.csv(text = txt, 
                                       stringsAsFactors = FALSE)
gauteng_weather_future_csv$Date <- as.Date(gauteng_weather_future_csv$Date)

# Convert the data frame into an xts object first (Date column as
# index), then into an idx_series anchored on the fitted model's own
# calendar so the resulting positions line up with res_weather$calendar.
gauteng_weather_future_xts <- xts::xts(
  gauteng_weather_future_csv[, -1],
  order.by = gauteng_weather_future_csv$Date
)
gauteng_weather_future_idx <- xts_to_idx(
  gauteng_weather_future_xts,
  start.pos = idx_to_pos(res_weather$calendar,
                         zoo::index(gauteng_weather_future_xts)[1])
)$series

# Supply CSV-based future xpred data to the model
# Keep the CSV-derived path in its own object and do NOT assign it to
# res_weather$xpred.new: res_weather$xpred.new already holds the
# full-precision `gauteng_weather_future` path set in section 3.1.2,
# and everything downstream (plot_log_forecast, plot_forecast,
# plot_compare_forecast) uses res_weather directly, so overwriting it
# here would silently replace full-precision weather with the
# 2-decimal-place CSV illustration for every later figure/table. This
# block exists only to demonstrate the CSV-reading pattern and to
# verify (via stopifnot) that the CSV scenario is consistent with the
# full-precision series - it does not feed into any later output.
res_weather_csv <- gauteng_weather_future_idx
stopifnot(
  nrow(res_weather_csv) == n.forecasts,
  isTRUE(all.equal(
    as.numeric(idx_values(res_weather_csv)[, "temperature_C"]),
    as.numeric(idx_values(gauteng_weather_future)[, "temperature_C"]),
    tolerance = 0.01
  ))
)
# res_weather$xpred.new is intentionally left untouched here (still the
# full-precision path from 3.1.2). If a user genuinely wants to forecast
# from the CSV-derived scenario instead, assign res_weather_csv to a
# *copy* of res_weather (or a dedicated results object), not to
# res_weather itself, so the full-precision run remains reproducible.

# Once xpred has been supplied, the fitted model will generate forecasts
# that are conditional on these external regressors.

# ----------------------------------------------------
# 3.1.4 Forecasts and Accuracy Plots (with Regressors)
# ----------------------------------------------------

# Forecast log growth from the weather-augmented model.
p <- tsgc::plot_log_forecast(
  res_weather,
  Y            = cumulative_cases,
  n.ahead      = n.forecasts,
  plt.start    = tail_pos_minus(res_weather$index, plt.length),
  title        = "Forecast of log growth rate of cases\n(with regressors: weather, oracle/realised)"
)
print(p)

# Forecast new cases conditional on the supplied future weather values.
p <- tsgc::plot_forecast(
  res_weather,
  n.ahead          = n.forecasts,
  confidence.level = CONF_LEVEL,
  plt.start        = tail_pos_minus(res_weather$index, plt.length),
  title            = "Forecast of new cases\nwith regressors (weather, oracle/realised), Gauteng",
  series.name      = "Cases"
)
print(p)

# Holdout accuracy: two weeks prior to end of sample
# Define a holdout estimation end 14 positions before est.end.1
# Reuse the same 14-day holdout design for the model with regressors.
est.end.holdout <- est.end.1 - n.forecasts

# Refit ONLY for holdout evaluation on the truncated window
# Refit the weather-augmented model on the truncated estimation sample.
model_q_xpred_holdout <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases,  
  xpred = gauteng_weather_est,
  q = q.default, 
  start = est.start.1,
  end   = est.end.holdout,
  calendar = gauteng_cal
)
# Estimate the truncated model before supplying holdout-period weather values.
res_qxpred_holdout <- tsgc::estimate(model_q_xpred_holdout)

# Extract weather values covering exactly the validation horizon.
gauteng_weather_holdout <- get_timeframe(
  gauteng_weather_idx,
  est.end.holdout + 1,
  est.end.holdout + n.forecasts
)

# Supply the holdout-period regressors required for conditional validation forecasts.
res_qxpred_holdout$xpred.new <- gauteng_weather_holdout

# Plot holdout accuracy for the model with weather regressors.
p <- tsgc::plot_holdout(
  res_qxpred_holdout,
  Y                = cumulative_cases,
  n.ahead          = n.forecasts,
  confidence.level = CONF_LEVEL,
  title            = "Accuracy: Forecast of new cases\nwith regressors (weather, oracle/realised)",
  series.name      = "Cases"
)
print(p)

# We can also compare different estimates with the actual trajectory
# Compare forecasts from the baseline fixed-q model and the weather-augmented model.
p <- tsgc::plot_compare_forecast(
  list(res_free, res_q, res_ar1, res_weather),
  actual = cumulative_cases
)
print(p)

# ============================
# 4. Reproduction Number (R_t)
# ============================

#' 
#' # 4. Reproduction Number (R_t)
#' 
#' Map Gompertz estimates to effective reproduction numbers using 
#' R_t = exp(g_y,t * gen_int), an assumed generation interval (4 days), 
#' reporting the most recent 7 daily estimates.
#' 
#' ## 4.1 Transform Gompertz estimates to R_t
#' 

## ---- 4.1 Transform Gompertz estimates to R_t ----
# Convert fitted Gompertz dynamics into an implied reproduction-number path using the assumed generation interval.
# Estimate the implied reproduction number from the fixed-q Gauteng model.
# estimate_r0() returns a dated data frame with fit/lower/upper columns
# over the last n.ahead (here ndays) positions; plotting is done explicitly
# with ggplot rather than a built-in show_plot argument.
r.t <- tsgc::estimate_r0(res_q, gen_int, ndays)

r.t.df <- data.frame(
  Date  = r.t$Date,
  Rt    = r.t$fit,
  lower = r.t$lower,
  upper = r.t$upper
)
names(r.t.df) <- c(
  "Date", "Rt",
  paste0("Rt_lower_", confidence_suffix(CONF_LEVEL)),
  paste0("Rt_upper_", confidence_suffix(CONF_LEVEL))
)

if (SAVE_TABLES) {
  # Save the R_t table for later use in reports or supplementary material.
  write.csv(r.t.df, row.names = FALSE, 
            file = file.path(tables_dir, "gauteng_gompertz_q005_rt.csv"))
  message("Saved gauteng_gompertz_q005_rt.csv")
}

# Build the R_t plot explicitly, since estimate_r0() does not take a
# show_plot/title argument.
fit <- lower <- upper <- NULL
p <- ggplot(data.frame(x = r.t.df$Date, fit = r.t.df$Rt,
                       lower = r.t.df[[3]], upper = r.t.df[[4]]),
            aes(x = x)) +
  geom_line(aes(y = fit, color = "Rt")) +
  geom_point(aes(y = fit), color = "red", size = 3) +
  geom_segment(aes(xend = x, yend = lower, y = fit), color = "blue") +
  geom_segment(aes(xend = x, yend = upper, y = fit), color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "68%  Interval"), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "solid", linewidth = 1.5, color = "black") +
  scale_x_date(date_breaks = "1 day", labels = scales::date_format("%d %b %y")) +
  labs(title = "Gauteng Reproduction numbers",
       x = "Date", y = expression(R[t])) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.2),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold")
  )
print(p)
save_plot(p, "gauteng_rt_gomp_q005_plot.png")

## ---- 4.2 Generation-interval sensitivity ----
# *** SENSITIVITY ANALYSIS: assumed generation interval (gen_int) ***
# R_t is derived from the fitted growth-rate path via a specified
# generation interval (gen_int = 4 days above); this value is an
# epidemiological assumption, not something estimated from this data,
# so the resulting R_t path is conditional on it. Rather than reporting
# a single R_t path without interrogating that assumption, sweep
# gen_int over a plausible range and compare the resulting R_t
# estimates on the same dates, so the sensitivity of the conclusion
# (e.g. whether R_t is above or below 1) to this assumption is visible
# rather than only asserted.
gen_int_grid <- c(3, 4, 5, 6, 7)

r.t.sensitivity <- lapply(gen_int_grid, function(g) {
  r.t.g <- tsgc::estimate_r0(res_q, g, ndays)
  data.frame(
    Date    = r.t.g$Date,
    gen_int = g,
    Rt      = r.t.g$fit
  )
})
r.t.sensitivity.df <- do.call(rbind, r.t.sensitivity)

if (SAVE_TABLES) {
  write.csv(r.t.sensitivity.df, row.names = FALSE,
            file = file.path(tables_dir, "gauteng_gompertz_q005_rt_gen_int_sensitivity.csv"))
  message("Saved gauteng_gompertz_q005_rt_gen_int_sensitivity.csv")
}

# Interpretation, not just the numbers: report the range of R_t implied
# across the assumed generation-interval grid on the most recent common
# date, and whether the above/below-1 conclusion is robust to it.
last_date <- max(r.t.sensitivity.df$Date)
rt_last <- r.t.sensitivity.df$Rt[r.t.sensitivity.df$Date == last_date]
message(
  "R_t sensitivity to gen_int on ", last_date, ": range [",
  paste(round(range(rt_last), 3), collapse = ", "), "] across gen_int in {",
  paste(gen_int_grid, collapse = ", "), "} days. ",
  if (all(rt_last > 1) || all(rt_last < 1)) {
    "The above/below-1 conclusion is unchanged across this grid."
  } else {
    "The above/below-1 conclusion CHANGES within this grid - R_t on this date is not robust to the generation-interval assumption."
  }
)

p_rt_sens <- ggplot2::ggplot(
  r.t.sensitivity.df,
  ggplot2::aes(x = Date, y = Rt, color = factor(gen_int))
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  ggplot2::labs(
    title = "Sensitivity of Gauteng R_t to the assumed generation interval",
    subtitle = paste0("gen_int swept over {", paste(gen_int_grid, collapse = ", "), "} days; point estimate only (no CI)"),
    x = "Date", y = expression(R[t]), color = "gen_int (days)"
  ) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    plot.title  = ggplot2::element_text(face = "bold")
  )
print(p_rt_sens)
save_plot(p_rt_sens, "gauteng_rt_gen_int_sensitivity.png")

# ========================================
# 5. Reinitialisation for Subsequent Waves
# ========================================

#' # 5. Reinitialisation for Subsequent Waves
#' 
#' Identify a trigger for reinitialising based on smoothed slope uncertainty; 
#' re-estimate with reinitialisation and assess forecasts vs baseline.
#' 
#' ## 5.1 Trigger Diagnostics
#' 
#' Compute smoothed slope and uncertainty; detect threshold crossings 
#' as reinitialisation triggers. 
#' The trigger is the first point where the estimated slope rises above 
#' its own 2-sigma upper confidence band after being below it 
#' in the previous period. 
#' If so, set the reinitialisation position at the latest sign-change prior 
#' to the first \(2\sigma\) crossing, marking the potential start 
#' of a new growth phase. 
#' 

## ---- 5.1 Reinitialisation Trigger Setup ----
# Estimate a longer Gauteng model and derive diagnostics used to identify possible reinitialisation dates.
# Extend the estimation window so the model can detect later-wave dynamics.
est.end.2 <- idx_to_pos(gauteng_cal, as.Date("2021-06-25"))

# Fit a baseline model over the longer window without reinitialisation.
model_rei_base <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, q = q.default,
  start = est.start.1, end = est.end.2,
  calendar = gauteng_cal
)
# Estimate the longer-window model and inspect the summary before deriving diagnostics.
res_rei_base <- tsgc::estimate(model_rei_base)
summary(res_rei_base)
tsgc::print_model_diagnostics(res_rei_base)

# KFS pieces from the fitted results object, extracted as idx_series so
# they carry their own integer positions.
smoothed.slope.full <- idx_series(res_rei_base$output$alphahat[, "slope"],
                                  start = res_rei_base$index[1])
# Use the smoothed-state covariance V (not the filtered Ptt/P), so the
# variance is matched to the smoothed alphahat estimate used above. This
# is a retrospective diagnostic, so the smoothed pair is the right one.
V.smoothed <- get_V(res_rei_base$output)
i.slope <- grep("slope", colnames(res_rei_base$output$alphahat))
smoothed.P.slope <- idx_series(V.smoothed[i.slope, i.slope, ],
                               start = res_rei_base$index[1])

# Combine slope estimates, uncertainty bands, and observed series into
# one diagnostic idx_series.
common_pos <- idx_positions(smoothed.slope.full)
d2 <- idx_series(
  cbind(
    smthd.slpe             = idx_values(smoothed.slope.full),
    plus.sd.smthd.slpe     = sqrt(idx_values(smoothed.P.slope)),
    plus.sd.smthd.slpe.1.5 = 1.5 * sqrt(idx_values(smoothed.P.slope)),
    plus.sd.smthd.slpe.2   = 2 * sqrt(idx_values(smoothed.P.slope))
  ),
  start = common_pos[1]
)

d2.mat <- as.matrix(idx_values(d2))
d2.df <- data.frame(
  Date = idx_to_date(gauteng_cal, idx_positions(d2)),
  smthd.slpe = d2.mat[, "smthd.slpe"],
  plus.sd.smthd.slpe = d2.mat[, "plus.sd.smthd.slpe"],
  plus.sd.smthd.slpe.1.5 = d2.mat[, "plus.sd.smthd.slpe.1.5"],
  plus.sd.smthd.slpe.2 = d2.mat[, "plus.sd.smthd.slpe.2"]
)
d2.df <- dplyr::filter(d2.df, Date >= as.Date("2020-10-06"))

# zt = smoothed slope minus its own two-standard-error threshold, i.e.
# the lower bound of the two-SE confidence interval around the slope.
# A trigger is a sign change in zt itself (zt > 0, lag(zt) <= 0), not a
# comparison of the lagged slope against the current period's threshold
# (which is not equivalent once the standard error changes over time).
d2.df$zt <- d2.df$smthd.slpe - d2.df$plus.sd.smthd.slpe.2

# Identify candidate dates where the two-SE lower bound crosses zero.
trigger.df <- d2.df %>%
  dplyr::mutate(prev_zt = dplyr::lag(zt)) %>%
  dplyr::filter(zt > 0 & prev_zt <= 0)

# Identify dates where the smoothed slope crosses or approaches zero, another restart diagnostic.
reinit_zero.df <- d2.df %>%
  dplyr::mutate(prev_smthd.slpe = dplyr::lag(smthd.slpe)) %>%
  dplyr::filter(Date < min(trigger.df$Date) & 
                  (smthd.slpe > 0 & prev_smthd.slpe < 0)) %>%
  dplyr::arrange(dplyr::desc(Date)) %>%
  dplyr::slice(1)

#' 
## ---- 5.1 Reinitialisation Trigger Plot ----
# Visualise the reinitialisation trigger criteria so the selected restart date is transparent.
# Build the diagnostic plot showing observed cases, slope behaviour, and candidate trigger dates.
p_trigger <-
  ggplot2::ggplot(data = d2.df[d2.df$Date > as.Date("2021-02-11"), ], 
                  ggplot2::aes(x = Date)) +
  ggplot2::geom_line(ggplot2::aes(y = smthd.slpe, color = "Smoothed slope"), 
                     linewidth = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe,
                                  color = "1 SE threshold"), 
                     linewidth = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe.1.5, 
                                  color = "1.5 SE threshold"), 
                     linewidth = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe.2, 
                                  color = "2 SE threshold"), 
                     linewidth = 0.5) +
  ggplot2::scale_y_continuous(n.breaks = 10) +
  ggplot2::geom_hline(yintercept = 0, linetype = "solid", 
                      color = "black", linewidth = 1) +
  ggplot2::geom_vline(data = trigger.df, 
                      ggplot2::aes(xintercept = Date), 
                      linewidth = 0.5, color = "black", linetype = "dashed") +
  ggplot2::geom_vline(data = reinit_zero.df, 
                      ggplot2::aes(xintercept = Date), 
                      linewidth = 1, color = "black") +
  ggplot2::labs(
    title = "Reinitialisation trigger diagnostic for Gauteng",
    subtitle = "Smoothed slope with 1, 1.5, and 2 standard-error lower thresholds (slope - k*se)",
    x = "Date",
    y = "Smoothed slope",
    caption = paste(
      "Thick vertical line: reset date. Dashed vertical line: 2-SE trigger date.",
      "Reset date is selected retrospectively from one episode/origin and is",
      "not validated as a general real-time reinitialisation rule."
    )
  ) +
  ggplot2::scale_x_date(date_breaks = "10 days") +
  ggplot2::scale_color_manual(
    name   = "Series", 
    values = c(
      "Smoothed slope" = "red",
      "1 SE threshold" = "blue",
      "1.5 SE threshold" = "green",
      "2 SE threshold" = "black"
    )
  ) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(
    legend.title = ggplot2::element_text(size = 9),
    legend.text  = ggplot2::element_text(size = 9),
    axis.text.x  = ggplot2::element_text(angle = 45, 
                                         hjust = 1, size = 8),
    plot.title   = ggplot2::element_text(face = "bold")
  )
print(p_trigger)
save_plot(p_trigger, "gauteng_cases_gomp_q005_reinit_trigger.png")

#' 
#' ## 5.2 Forecasts & Accuracy (post-reinitialisation)
#' 
#' Re-estimate with a chosen reinit position; produce forecasts 
#' and compare to the no-reinit baseline.
#' 

## ---- 5.2 Reinitialisation Estimation & Forecasts ----
# Refit the Gauteng model with a reinitialisation position and compare forecasts with and without reinitialisation.
# Set reinit position (could also take from trigger.df/reinit_zero.df)
# This chosen date starts the model state afresh for the later wave.
reinit.pos <- idx_to_pos(gauteng_cal, as.Date("2021-04-21"))

# Fit the dynamic Gompertz model with reinitialisation activated at the selected position.
model_reinit <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, q = q.default,
  start = est.start.1, end = est.end.2,
  reinit.idx = reinit.pos,
  calendar = gauteng_cal
)
# Estimate the reinitialised model and inspect its summary.
res_reinit <- tsgc::estimate(model_reinit)
summary(res_reinit)
tsgc::print_model_diagnostics(res_reinit)

# Forecasts after reinitialisation
# Forecast log growth after allowing the model to restart at the reinitialisation position.
p <- tsgc::plot_log_forecast(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start = tail_pos_minus(res_reinit$index, plt.length),
  title = "Forecast of log growth rate of cases\nafter reinitialisation"
); print(p)

# Forecast new cases from the reinitialised model.
p <- tsgc::plot_forecast(
  res_reinit, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start = tail_pos_minus(res_reinit$index, plt.length),
  title = "Forecast of new cases\nafter reinitialisation", 
  series.name = "Cases"
); print(p)

# Evaluate holdout accuracy for the reinitialised model.
p <- tsgc::plot_holdout(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL, 
  title = "Accuracy: Forecast of new cases\nwith reinitialisation", 
  series.name = "Cases"
); print(p)

# Compare holdout with/without reinitialisation (baseline = res_rei_base)
# Produce the comparable holdout plot for the non-reinitialised longer-window model.
tsgc::plot_holdout(
  res_rei_base, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL, 
  title = "Accuracy: Forecast of new cases\nwithout reinitialisation",
  series.name = "Cases"
)
# Compare the reinitialised and non-reinitialised forecasts directly.
tsgc::plot_compare_forecast(
  list(res_rei_base, res_reinit), 
  actual = cumulative_cases
)

# ==========================
# 6. Leading Indicator Model
# ==========================

#' 
#' # 6. Leading Indicator Model: England (Daily)
#' 
#' Model how cases (lead) anticipate hospitalisations (target) at lag L, 
#' with weekly seasonality; produce forecasts and accuracy metrics. 
#' 
#' ## 6.1 Leading Indicator: Baseline Model
#' 
#' Estimate a parsimonious leading-indicator model 
#' and evaluate short-horizon performance.
#' 

# ---- 6.1 Baseline: England ----
# Fit the England leading-indicator example, where one series helps forecast the target hospital-admissions series.
# Load and convert the England data (cases, hospitalisations) to idx_series/idx_calendar form.
data(england, package = "tsgc")
conv <- xts_to_idx(england[, 1:2])
eng <- conv$series
eng_cal <- conv$calendar

# Quick plot
# Create a quick leading-indicator model to visualise the lead/target relationship with n.lag = 4.
mod2 <- tsgc::SSModelLeadingIndicator(eng, n.lag = 4, calendar = eng_cal)
p <- plot(
  mod2, title = "Daily COVID cases and Hospitalisations\n(England)",
  series.name.lead = "Cases", series.name.target = "Hospitalisations", 
  take.log = TRUE
)
print(p)
save_plot(p, "eng_hosp_lead_cases.png")

# Estimation window
# Define the England estimation window, plotting length, lag, and forecast horizon.
est.start.eng <- idx_to_pos(eng_cal, "2021-04-30")
est.end.eng   <- idx_to_pos(eng_cal, "2021-07-24")
plt.len.eng   <- 14
n.lag         <- 4
n.forecasts   <- 7

# Define and estimate
# Specify the England leading-indicator model for hospital admissions.
out_eng <- tsgc::SSModelLeadingIndicator(
  Y = eng, n.lag = n.lag, q = NULL, LeadIndCol = 1, sea.period = 7,
  start = est.start.eng, end = est.end.eng, calendar = eng_cal
)
# Estimate the England leading-indicator model and inspect the fitted summary.
res_eng <- tsgc::estimate(out_eng)
summary(res_eng)
tsgc::print_model_diagnostics(res_eng)

# Forecasts
# Forecast the target-series log growth rate.
p <- tsgc::plot_log_forecast(
  res_eng, Y = eng, n.ahead = n.forecasts, 
  plt.start = est.end.eng - plt.len.eng,
  title = "Forecast of log growth rate of hospital admissions\n(England)"
); print(p)

# Forecast hospital admissions in the original scale.
p <- tsgc::plot_forecast(
  res_eng, n.ahead = n.forecasts, 
  plt.start = est.end.eng - plt.len.eng,
  series.name = "Hospital admissions", 
  title = "Forecast of hospital admissions\n(England)"
); print(p)

# Evaluate the England forecast against a holdout segment.
p <- tsgc::plot_holdout(
  res_eng, Y = eng, n.ahead = n.forecasts,
  series.name = "Hospital admissions", 
  title = "Accuracy: Forecast of hospital admissions\n(England)"
); print(p)

if (SAVE_TABLES) {
  # Export the England leading-indicator forecast and filtered-state outputs as CSV files.
  write_results_clear(
    res = res_eng,
    res.dir = tables_dir,
    n.ahead = n.forecasts,
    model_slug = "england_leading_indicator",
    target_slug = "hospital_admissions",
    confidence.level = CONF_LEVEL
  )
  write_csv_manifest(tables_dir)
  message("Saved clear CSV results for: england_leading_indicator")
}

#' 
#' ## 6.2 Leading Indicator: Model with Regressors
#' 
#' Augment with exogenous weather regressors for both lead and target, 
#' then forecast and evaluate.
#' 

## ---- 6.2 England With Regressors - xpred, eval=TRUE ----
# Extend the England leading-indicator model with weather regressors for both lead and target series.
# Load and convert the weather regressors, anchored on the England calendar.
data(england_weather_2021, package = "tsgc")
conv <- xts_to_idx(
  england_weather_2021[, 1:4],
  start.pos = idx_to_pos(eng_cal, zoo::index(england_weather_2021)[1])
)
england_weather_idx <- conv$series
xpred_lead <- xpred_targ <- england_weather_idx

# Fit the England leading-indicator model with xpred regressors.
mod_eng_x <- tsgc::SSModelLeadingIndicator(
  eng, n.lag = 4, xpred_lead = xpred_lead, xpred_targ = xpred_targ,
  start = est.start.eng, end = est.end.eng, calendar = eng_cal
)
# Estimate the weather-augmented England model and inspect the summary.
res_eng_x <- tsgc::estimate(mod_eng_x)
summary(res_eng_x)
tsgc::print_model_diagnostics(res_eng_x)

# Supply future regressors for the lead and target equations. 
res_eng_x$xpred_lead.new <- england_weather_idx
res_eng_x$xpred_targ.new <- england_weather_idx

# Forecast log growth from the weather-augmented leading-indicator model.
p <- tsgc::plot_log_forecast(
  res_eng_x, Y = eng, n.ahead = n.forecasts, 
  plt.start = est.end.eng - plt.len.eng,
  title = "Forecast of log growth rate of hospital admissions\nwith regressors (weather, oracle/realised), England"
); print(p)

# Forecast hospital admissions with weather regressors included.
p <- tsgc::plot_forecast(
  res_eng_x, n.ahead = n.forecasts, 
  plt.start = est.end.eng - plt.len.eng,
  title = "Forecast of hospital admissions\nwith regressors (weather, oracle/realised), England", 
  series.name = "Hospital admissions"
); print(p)

# Evaluate holdout accuracy for the weather-augmented England model.
p <- tsgc::plot_holdout(
  res_eng_x, Y = eng, n.ahead = n.forecasts,
  title = "Accuracy: Forecast of hospital admissions\nwith regressors (weather, oracle/realised), England", 
  series.name = "Hospital admissions"
); print(p)

# ===================================================
# 7. Comparing Leading Indicator and Gompertz Models
# ===================================================

#' 
#' # 7. Comparing Leading Indicator and Gompertz Models: UK & Italy
#' 
#' Compare UK forecasts from a UK-only Gompertz model 
#' vs a leading-indicator model with Italy as the lead across two windows. 
#' Both models are estimated on the same window with horizon 14 days 
#' and use the same confidence level.
#' 
#' ## 7.1 Case 1: First Peak (UK vis-a-vis Italy)
#' 

# ---- 7.1 UK vis-a-vis Italy ----
# Compare a UK-only Gompertz forecast with an Italy-to-UK leading-indicator forecast during the first peak window.

# By default, column 1 is treated as the leading indicator and 
# column 2 as the target.

# Load and convert the UK-Italy data to idx_series/idx_calendar form.
data(ukitaly, package = "tsgc")
conv <- xts_to_idx(ukitaly)
ukitaly_idx <- conv$series
ukitaly_cal <- conv$calendar

# Create a quick UK-Italy leading-indicator object to inspect the lead/target timing.
ukit <- tsgc::SSModelLeadingIndicator(ukitaly_idx, n.lag = 4, calendar = ukitaly_cal)
p <- plot(
  ukit, title = "Daily COVID cases in UK and Italy",
  series.name.lead   = "Italy",
  series.name.target = "UK", 
  take.log = FALSE
)
print(p)
save_plot(p, "ukit_cases_lead_cases.png")

# 7.1 Case 1: First peak window
# Set the forecast horizon, plotting window, and confidence level for the UK-Italy comparisons.
n.forecasts <- 14
plt.length  <- 30
CONF        <- CONF_LEVEL

# Case 1 uses the first-peak estimation window ending 1 April 2020.
est.start <- idx_to_pos(ukitaly_cal, "2020-02-25")
est.end   <- idx_to_pos(ukitaly_cal, "2020-04-01")
Yuk       <- idx_series(idx_values(ukitaly_idx)[, "UK"], start = ukitaly_idx$start)

# Estimate a UK-only dynamic Gompertz benchmark.
res_uk_gomp1 <- tsgc::estimate(
  tsgc::SSModelDynamicGompertz(
    Y = Yuk, q = q.default,
    start = est.start, 
    end   = est.end,
    calendar = ukitaly_cal
  )
)

# Forecast UK daily cases from the Gompertz benchmark.
p <- tsgc::plot_forecast(
  res_uk_gomp1, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz)",
  plt.start = tail_pos_minus(res_uk_gomp1$index, plt.length),
  series.name    = "UK cases"
); print(p)

# Evaluate holdout accuracy for the UK Gompertz benchmark.
p <- tsgc::plot_holdout(
  res_uk_gomp1, Y = Yuk, n.ahead = n.forecasts, 
  confidence.level = CONF,
  title      = "Accuracy: Forecast of daily COVID cases\nUK (Gompertz)", 
  series.name = "UK cases"
); print(p)

# Leading indicator model, Case 1
# Use a 14-day Italy-to-UK lag for the leading-indicator comparison.
# NOTE: this fixed value is illustrative only; it is not the output of
# the systematic lag comparison performed later (Section 7.2).
n.lag <- 14
# Estimate the Italy-to-UK leading-indicator model for the same window.
res_uk_lead1 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(
    Y = ukitaly_idx, 
    n.lag = n.lag, sea.period = 7,
    start = est.start, 
    end   = est.end,
    calendar = ukitaly_cal
  )
)

# Forecast UK daily cases from the leading-indicator model.
p <- tsgc::plot_forecast(
  res_uk_lead1, n.ahead = n.forecasts,
  title = "Leading indicator forecast\ndaily COVID cases in UK",
  plt.start = est.end - 30, series.name = "UK cases"
); print(p)

# Evaluate holdout accuracy for the leading-indicator forecast.
p <- tsgc::plot_holdout(
  res_uk_lead1, Y = ukitaly_idx, n.ahead = n.forecasts,
  title      = "Accuracy: Leading indicator forecast\ndaily COVID cases in UK", 
  series.name = "UK cases"
); print(p)

# Plot the Gompertz and leading-indicator forecasts together for visual comparison.
# actual must be the single UK column (Yuk), not the full two-column
# ukitaly_idx: plot_compare_forecast()'s actual.diff is built assuming a
# univariate series, and a two-column series here causes an rbind()
# column-count mismatch when combining forecasts with the actual values.
tsgc::plot_compare_forecast(
  list(res_uk_gomp1, res_uk_lead1), 
  actual = Yuk
)

# -----------------------------
# 7.2 Cross-Validation (UK/IT)
# -----------------------------

#' 
#' ## 7.2 Cross-Validation: UK Gompertz vs Leading-Indicator Models
#' 
#' Compare out-of-sample performance via rolling cross-validation across
#' a set of Gompertz and leading-indicator specifications.
#' 

## ---- 7.2 Cross Validation ----
# Run cross-validation over baseline Gompertz models and leading-indicator lag choices to compare forecast performance systematically.
#
# Methodological note (see accompanying critical assessment): the
# original design compared 23 candidates using 5 closely spaced
# forecast origins (gap = 2) with a 14-day horizon, so consecutive
# origins' evaluation windows overlapped heavily, and the same folds
# were used both to select the winning lag and to report its accuracy.
#
# cross_val() steps origins FORWARD from est.end by gap each iteration
# (model$end <- est.end + (k-1)*gap), not backward. This block
# therefore lays out two disjoint, forward-stepping blocks inside the
# fixed window [est.start, est.end.cv]: an earlier SELECTION block
# (used only to pick the lag) and a later, non-overlapping REPORTING
# block (used only to score the selected model), with the REPORTING
# block's final forecast horizon landing exactly at est.end.cv so
# nothing runs past the available UK-Italy sample, and the SELECTION
# block's last forecast horizon finishing strictly before (not merely
# up to) the REPORTING origin - verified below by an explicit
# stopifnot(), not asserted by comment alone. A 7-day horizon (rather
# than 14) is used so both blocks fit inside the available ~50-day
# window while still leaving a reasonable estimation sample before the
# first SELECTION origin.
n.ahead.cv    <- 7
est.end.cv    <- idx_to_pos(ukitaly_cal, "2020-04-15")  # Case 2's extended window end
gap.cv        <- n.ahead.cv                              # non-overlapping horizons
n.select.cv   <- 2
n.report.cv   <- 1

# REPORTING block: est.end chosen so its last origin's n.ahead.cv-day
# forecast horizon lands exactly at est.end.cv (no overrun past the
# available data).
report.end.cv <- est.end.cv - (n.report.cv - 1) * gap.cv - n.ahead.cv
# SELECTION block: est.end chosen so its origins and their forecast
# horizons finish strictly before the REPORTING block's first origin
# (disjoint, non-overlapping folds).
#
# The last SELECTION origin is select.end.cv + (n.select.cv-1)*gap.cv,
# with its forecast horizon reaching select.end.cv + (n.select.cv-1)*
# gap.cv + n.ahead.cv. The REPORTING block's first (only) origin is
# report.end.cv itself. Requiring these strictly non-overlapping means:
#   select.end.cv + (n.select.cv - 1)*gap.cv + n.ahead.cv < report.end.cv
# i.e.
#   select.end.cv < report.end.cv - (n.select.cv - 1)*gap.cv - n.ahead.cv
# Subtracting n.select.cv * gap.cv (rather than a single gap.cv) leaves
# a full extra gap between the last SELECTION horizon and the REPORTING
# origin, giving strict separation for any n.select.cv/n.report.cv.
select.end.cv <- report.end.cv - n.select.cv * gap.cv - n.ahead.cv
stopifnot(select.end.cv > est.start)

# Explicit disjointness check (not just an arithmetic comment): the
# last SELECTION fold's forecast horizon must end strictly before the
# REPORTING block's first origin.
last_selection_origin      <- select.end.cv + (n.select.cv - 1) * gap.cv
last_selection_horizon_end <- last_selection_origin + n.ahead.cv
stopifnot(last_selection_horizon_end < report.end.cv)

cv_models <- list()

# Naive benchmarks
cv_models[["Naive_last_value"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, q = 0, start = est.start, end = est.end.cv, calendar = ukitaly_cal
)
cv_models[["RW_growth"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, q = q.default, sea.period = 0, start = est.start,
  end = est.end.cv, calendar = ukitaly_cal
)

# Model 1: Vanilla Gompertz
cv_models[["Vanilla_q"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, q = q.default, start = est.start, 
  end = est.end.cv, calendar = ukitaly_cal
)

# Model 2: Vanilla Gompertz with AR(1)
cv_models[["Vanilla_ar1"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, start = est.start, 
  end = est.end.cv, ar1 = TRUE, calendar = ukitaly_cal
)

# Leading-indicator lag candidates
# NOTE: the 20 lag is the minimum for two rows to be available.
# This leaves a very thin estimation sample. One could opt to
# use a later selection block or relax the full 7-day buffer 
# before reporting.
for (i in 1:20) {
  cv_models[[paste0("Lag", i)]] <- tsgc::SSModelLeadingIndicator(
    Y = ukitaly_idx, start = est.start, 
    end = est.end.cv, n.lag = i, calendar = ukitaly_cal
  )
}

# SELECTION run: earlier, disjoint origins only, used to choose the lag.
cv_selection <- tsgc::cross_val(
  Y           = ukitaly_idx, 
  model_list  = cv_models, 
  est.end     = select.end.cv,
  criterion   = "smape",
  n.ahead     = n.ahead.cv,
  n.estimate  = n.select.cv, 
  gap         = gap.cv
)
print(cv_selection)

# Pick the best-performing Lag* candidate on the SELECTION folds only.
lag_rows <- cv_selection[grepl("^Lag", cv_selection$Model), ]
lag_rows$mean_smape <- rowMeans(lag_rows[, -1, drop = FALSE])
selected_lag_name <- lag_rows$Model[which.min(lag_rows$mean_smape)]
message("Lag selected on SELECTION folds: ", selected_lag_name)

# REPORTING run: later, disjoint origin(s), used only to score the
# selected lag against the benchmarks and other candidates. This is
# the number that should be reported as the model's out-of-sample
# accuracy, not the SELECTION-fold number used to pick the lag.
cv_report_models <- cv_models[c(
  "Naive_last_value", "RW_growth", "Vanilla_q", "Vanilla_ar1", selected_lag_name
)]
cv_reporting <- tsgc::cross_val(
  Y           = ukitaly_idx, 
  model_list  = cv_report_models, 
  est.end     = report.end.cv,
  criterion   = "smape",
  n.ahead     = n.ahead.cv,
  n.estimate  = n.report.cv, 
  gap         = gap.cv
)
print(cv_reporting)

if (SAVE_TABLES) {
  write.csv(cv_selection, file.path(tables_dir, "cv_lag_selection_smape.csv"), row.names = FALSE)
  write.csv(cv_reporting, file.path(tables_dir, "cv_reporting_smape.csv"), row.names = FALSE)
}

# -------------------------------
# 7.3 Case 2: Extended Window
# -------------------------------

#' 
#' ## 7.3 Case 2: Extended Window (UK vis-a-vis Italy)
#' 
#' Re-estimate both models on an extended sample window and 
#' compare forecasts and accuracy.
#' 

# Case 2: Extended window
# Case 2 extends the estimation window to 15 April 2020 to test sensitivity to a longer sample.
est.start <- idx_to_pos(ukitaly_cal, "2020-02-25")
est.end   <- idx_to_pos(ukitaly_cal, "2020-04-15")

# Estimate the extended-window UK Gompertz model.
res_uk_gomp2 <- tsgc::estimate(
  tsgc::SSModelDynamicGompertz(
    Y = Yuk, q = q.default,
    start = est.start, end = est.end, calendar = ukitaly_cal
  )
)

# Forecast UK cases from the extended-window Gompertz model.
p <- tsgc::plot_forecast(
  res_uk_gomp2, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz, extended)",
  plt.start = tail_pos_minus(res_uk_gomp2$index, plt.length),
  series.name    = "UK cases"
); print(p)

# Evaluate holdout accuracy for the extended-window Gompertz model.
p <- tsgc::plot_holdout(
  res_uk_gomp2, Y = Yuk, n.ahead = n.forecasts, 
  confidence.level = CONF,
  title      = "Accuracy: Forecast of daily COVID cases\nUK (Gompertz, extended)", 
  series.name = "UK cases"
); print(p)

# Estimate the extended-window Italy-to-UK leading-indicator model.
# NOTE: n.lag = 14 here is illustrative and precedes the systematic lag
# comparison in Section 7.2, which selects its lag from data (see
# `selected_lag_name`, generally not 14) using disjoint
# selection/reporting folds. This fixed value is not justified by that
# comparison; it demonstrates the API only.
res_uk_lead2 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(
    Y = ukitaly_idx, n.lag = 14,
    start = est.start, end = est.end, calendar = ukitaly_cal
  )
)

# Forecast UK cases from the extended-window leading-indicator model.
p <- tsgc::plot_forecast(
  res_uk_lead2, n.ahead = n.forecasts,
  title = "Forecast of daily COVID cases\nUK (Leading indicator model, extended)",
  plt.start = est.end - plt.length, 
  series.name = "UK cases"
); print(p)

# Evaluate holdout accuracy for the extended-window leading-indicator model.
p <- tsgc::plot_holdout(
  res_uk_lead2, Y = ukitaly_idx, n.ahead = n.forecasts,
  title      = "Accuracy: Forecast of daily COVID cases\nUK (Leading indicator model, extended)", 
  series.name = "UK cases"
); print(p)

# Compare the two extended-window forecasts directly.
tsgc::plot_compare_forecast(
  list(res_uk_gomp2, res_uk_lead2), 
  actual = idx_series(idx_values(ukitaly_idx)[, "UK"], start = ukitaly_idx$start)
)

# ==========================
# 8. Other Data Frequencies
# ==========================

#' 
#' # 8. Extensions to Other Data Frequencies 
#' 
#' Demonstrate quarterly, monthly, and annual use-cases 
#' with appropriate seasonal settings and lead-lag structures.
#' 
#' ## 8.1 Quarterly: Wii
#' 
#' Sales series for the Nintendo Wii console. 
#' Gompertz model applied to a non-epidemic diffusion process observed quarterly.
#' 
#' Fit a quarterly Gompertz to Wii cumulative sales; forecast and evaluate.
#' 

## ---- 8.1 Quarterly: Wii ----
# Demonstrate that the same modelling interface works for quarterly sales data, not only daily epidemiological data.
# Load Nintendo quarterly sales data and build a quarterly idx_series/idx_calendar (idx_calendar directly, since xts_to_idx() assumes a daily step).
data(nintendo_sales, package = "tsgc")

nintendo_idx <- idx_series(zoo::coredata(nintendo_sales), start = 1L)
nintendo_cal <- idx_calendar(
  anchor = as.Date(zoo::index(nintendo_sales)[1]),
  anchor_pos = 1L,
  amount = 1, unit = "quarters",
  posixct = TRUE
)
wii <- idx_series(idx_values(nintendo_idx)[, "wii"], start = nintendo_idx$start)

# Note: The tsgc function requires the input cumulative series to be strictly 
# increasing in time. If the cumulative values exhibit plateaus—as in 
# the case of the Wii series—it is necessary to add small increments 
# to eliminate flat segments and allow model estimation. Note that a 
# strictly increasing cumulative series also implies that the underlying 
# (non-cumulative) series must be strictly positive.

# Model estimated for a strictly increasing segment
# Use a four-quarter forecast horizon and quarterly positions for the quarterly model.
n.forecasts <- 4
est.start.q <- idx_to_pos(nintendo_cal, "2006-10-01")
est.end.q   <- idx_to_pos(nintendo_cal, "2010-07-01")

# Fit a quarterly dynamic Gompertz model to Wii sales.
mod_wii <- tsgc::SSModelDynamicGompertz(  
  Y = wii, sea.period = 4, start = est.start.q, 
  end = est.end.q, calendar = nintendo_cal
)
# Estimate the quarterly Wii model and inspect the summary.
res_wii <- tsgc::estimate(mod_wii)
summary(res_wii)
tsgc::print_model_diagnostics(res_wii)

# Cases with MA overlay
p <- plot(
  mod_wii, title = "Quarterly Wii console sales", 
  series.name = "Sales (million units)", MA_period = 4
)
print(p)
save_plot(p, "wii_sales_gomp_cases.png")

# Forecast Wii log growth at quarterly frequency.
p <- tsgc::plot_log_forecast(
  res_wii, Y = wii, n.ahead = n.forecasts, 
  title = "Forecast of log growth rate of Wii sales"
)
print(p)
save_plot(p, "wii_sales_gomp_loggr_fcst.png")

# Forecast quarterly Wii sales in the original scale.
p <- tsgc::plot_forecast(
  res_wii, n.ahead = n.forecasts,
  title = "Forecast of new Wii sales",
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "wii_sales_gomp_fcst.png")

# Evaluate quarterly holdout accuracy for Wii sales.
# n=4 holdout observations (quarterly): treat coverage/accuracy here as
# illustrative of the API, not as a calibration claim (see #23 caveat).
p <- tsgc::plot_holdout(
  res_wii, Y = wii, n.ahead = n.forecasts, 
  title = paste0("Accuracy: Forecast of new Wii sales\n",
                 "(illustrative only - n=", n.forecasts, " holdout observations)"),
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "wii_sales_gomp_holdout.png")

#' ## 8.2 Leading Indicator Model with Quarterly Data: Wii to Switch
#' 
#'  Nintendo Wii was launched in 2006, and the Nintendo Switch in 2017.
#' 
#' Model quarterly lead-lag between Wii and Switch; forecast and evaluate.
#' 

## ----  8.2 Quarterly: Wii to Switch (Lead) ----
# Use Wii sales as a quarterly leading indicator for Switch sales and forecast the target series.
# Set the quarterly leading-indicator horizon and estimation window for Switch sales.
n.forecasts   <- 8
est.start.q2  <- idx_to_pos(nintendo_cal, "2017-01-01")
est.end.q2    <- idx_to_pos(nintendo_cal, "2019-10-01")
# Compute the lead-lag distance between Wii and Switch launches in quarters.
n.lag.q       <- idx_to_pos(nintendo_cal, "2017-01-01") - idx_to_pos(nintendo_cal, "2006-10-01")

# Column 1 (Wii) is treated as the lead; column 2 (Switch) as the target.
# Select the lead and target quarterly sales series.
y_q <- idx_series(idx_values(nintendo_idx)[, c("wii", "switch_all")], start = nintendo_idx$start)
# Fit the quarterly Wii-to-Switch leading-indicator model.
mod_switch <- tsgc::SSModelLeadingIndicator(
  Y = y_q, sea.period = 4, n.lag = n.lag.q, 
  start = est.start.q2, end = est.end.q2, calendar = nintendo_cal
)
# Estimate the quarterly leading-indicator model and inspect the summary.
res_switch <- tsgc::estimate(mod_switch)
summary(res_switch)
tsgc::print_model_diagnostics(res_switch)

# Forecast Switch log growth at quarterly frequency.
p <- tsgc::plot_log_forecast(
  res_switch, Y = y_q, n.ahead = n.forecasts, 
  title = "Forecast of log growth rate of Switch sales"
)
print(p)
save_plot(p, "switch_sales_lead_loggr_fcst.png")

# Forecast quarterly Switch sales.
p <- tsgc::plot_forecast(
  res_switch, n.ahead = n.forecasts, 
  title = "Forecast of new Switch sales",
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "switch_sales_lead_fcst.png")

# Evaluate holdout accuracy for the Switch forecast.
# n=8 holdout observations (quarterly): illustrative of the API, not a
# calibration claim (see #23 caveat).
p <- tsgc::plot_holdout(
  res_switch, Y = y_q, n.ahead = n.forecasts, 
  title = paste0("Accuracy: Forecast of new Switch sales\n",
                 "(illustrative only - n=", n.forecasts, " holdout observations)"),
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "switch_sales_lead_holdout.png")

#' 
#' ## 8.3 Monthly: Plus500
#' 
#' Plus500 is an online retail trading platform (fintech).
#' 
#' Fit a monthly Gompertz model to Plus500 app downloads; forecast and evaluate.
#' 

## ---- 8.3 Monthly-etrading: Plus500 ----
# Apply the dynamic Gompertz model to monthly app-download data.
# Load monthly e-trading app downloads and build a monthly idx_series/idx_calendar.
data(etrading_apps, package = "tsgc")

etrading_idx <- idx_series(zoo::coredata(etrading_apps), start = 1L)
etrading_cal <- idx_calendar(
  anchor = as.Date(zoo::index(etrading_apps)[1]),
  anchor_pos = 1L,
  amount = 1, unit = "months",
  posixct = TRUE
)
Plus500 <- idx_series(idx_values(etrading_idx)[, 1], start = etrading_idx$start)

# Use a four-month horizon and monthly positions for the monthly model.
n.forecasts <- 4
est.start.m <- idx_to_pos(etrading_cal, "2016-01-01")
est.end.m   <- idx_to_pos(etrading_cal, "2021-01-01")

# Fit a monthly dynamic Gompertz model to Plus500 downloads.
mod_500 <- tsgc::SSModelDynamicGompertz(
  Y = Plus500, sea.period = 12, start = est.start.m, 
  end = est.end.m, calendar = etrading_cal
)
# Estimate the monthly Plus500 model and inspect the summary.
res_500 <- tsgc::estimate(mod_500)
summary(res_500)
tsgc::print_model_diagnostics(res_500)

p <- plot(
  mod_500, title = "Plus500 monthly downloads in France", 
  series.name = "Monthly downloads", MA_period = 4
)
print(p)
save_plot(p, "plus500_downloads_gomp_cases.png")

# Forecast Plus500 log growth at monthly frequency.
p <- tsgc::plot_log_forecast(
  res_500, Y = Plus500, n.ahead = n.forecasts, 
  title = "Forecast of log growth rate of Plus500 downloads"
)
print(p)
save_plot(p, "plus500_downloads_gomp_loggr_fcst.png")

# Forecast monthly Plus500 downloads in the original scale.
p <- tsgc::plot_forecast(
  res_500, n.ahead = n.forecasts, 
  title = "Forecast of new Plus500 downloads",
  series.name = "Downloads"
)
print(p)
save_plot(p, "plus500_downloads_gomp_fcst.png")

# Evaluate holdout accuracy for the Plus500 forecast.
# n=4 holdout observations (monthly): illustrative of the API, not a
# calibration claim (see #23 caveat).
p <- tsgc::plot_holdout(
  res_500, Y = Plus500, n.ahead = n.forecasts, 
  title = paste0("Accuracy: Forecast of new Plus500 downloads\n",
                 "(illustrative only - n=", n.forecasts, " holdout observations)"),
  series.name = "Downloads"
)
print(p)
save_plot(p, "plus500_downloads_gomp_holdout.png")

#' 
#' ## 8.4 Leading Indicator Model with Monthly Data: DEGIRO to AvaTrade
#' 
#' This example uses two online retail trading apps.  
#' DEGIRO is a major European low-cost brokerage platform. 
#' AvaTrade is a global CFD/FX trading platform. 
#' 
#' Monthly leading-indicator model with DEGIRO as the lead for AvaTrade; 
#' forecast and evaluate.
#' 

## ---- 8.4 Monthly-leading: DEGIRO to AvaTrade (Lead) ----
# Use monthly DEGIRO downloads as a leading indicator for AvaTrade downloads.
# Set the monthly leading-indicator horizon, positions, and lag for AvaTrade.
n.forecasts  <- 4
est.start.m2 <- idx_to_pos(etrading_cal, "2017-07-01")
est.end.m2   <- idx_to_pos(etrading_cal, "2021-02-01")
# Compute the DEGIRO-to-AvaTrade lag in months.
n.lag.m      <- idx_to_pos(etrading_cal, "2017-07-01") - idx_to_pos(etrading_cal, "2017-01-01")

# Select the monthly lead and target app-download series.
y_m <- idx_series(idx_values(etrading_idx)[, c("DEGIRO", "AvaTrade")], start = etrading_idx$start)
# Fit the monthly DEGIRO-to-AvaTrade leading-indicator model.
mod_500_lead <- tsgc::SSModelLeadingIndicator(
  Y = y_m, sea.period = 12, n.lag = n.lag.m, 
  start = est.start.m2, end = est.end.m2, calendar = etrading_cal
)
# Estimate the monthly leading-indicator model and inspect the summary.
res_500_lead <- tsgc::estimate(mod_500_lead)
summary(res_500_lead)
tsgc::print_model_diagnostics(res_500_lead)

# Forecast AvaTrade log growth at monthly frequency.
p <- tsgc::plot_log_forecast(
  res_500_lead, Y = y_m, n.ahead = n.forecasts, 
  title = "Forecast of log growth rate of AvaTrade downloads"
)
print(p)
save_plot(p, "avatrade_downloads_lead_loggr_fcst.png")

# Forecast monthly AvaTrade downloads.
p <- tsgc::plot_forecast(
  res_500_lead, n.ahead = n.forecasts, 
  title = "Forecast of new AvaTrade downloads", 
  series.name = "Downloads"
)
print(p)
save_plot(p, "avatrade_downloads_lead_fcst.png")

# Evaluate holdout accuracy for the AvaTrade forecast.
# n=4 holdout observations (monthly): illustrative of the API, not a
# calibration claim (see #23 caveat).
p <- tsgc::plot_holdout(
  res_500_lead, Y = y_m, n.ahead = n.forecasts, 
  title = paste0("Accuracy: Forecast of new AvaTrade downloads\n",
                 "(illustrative only - n=", n.forecasts, " holdout observations)"), 
  series.name = "Downloads"
)
print(p)
save_plot(p, "avatrade_downloads_lead_holdout.png")

#' 
#' ## 8.5 Annual: 3DS
#' 
#' 3DS is Nintendo's handheld console (released 2011).
#' Here we subsample quarterly cumulative global sales to an
#' annual-frequency series, taking every 4th quarter starting from the
#' first observation in the source series, to demonstrate annual
#' Gompertz modelling. These are same-quarter year-over-year
#' cumulative-sales differences (whatever quarter the source series
#' happens to start on), not calendar-year sales.
#' 
#' Annual modelling with `sea.period = 0`; 
#' forecast and evaluate 3DS annual series.
#' 

## ---- 8.5 Annual: 3DS ----
# Aggregate quarterly Nintendo data to annual frequency and fit an annual Gompertz model.
# Use a two-year forecast horizon and annual positions for the annual model.
n.forecasts <- 2

# Subsample quarterly Nintendo sales to annual-frequency observations
# (every 4th quarter, i.e. the same calendar quarter each year). We
# preserve the true selected dates rather than relabelling them onto a
# 1 January calendar, so the resulting series is accurately described
# as annual-frequency same-quarter cumulative sales, not calendar-year
# sales. Rather than hardcoding an assumed starting month (which
# depends on which quarter the source series happens to begin on and
# was simply wrong here - it does not fall in Jul/Aug/Sep), we verify
# directly against the series' own first date: every selected date
# must fall in the same calendar month as the first observation, since
# a 4-quarter step always returns to the same quarter.
# zoo::index(nintendo_sales) is a `yearqtr` object, not a `Date`. Every
# downstream use here (idx_calendar()'s anchor, and the format()/
# idx_to_pos() calls further below) needs a genuine Date - passing a
# yearqtr straight through as an idx_calendar anchor with unit = "years"
# errors ("a yearqtr cal$anchor requires cal$unit = 'quarters'"), and
# format()ing a yearqtr with "%Y-%m-%d" does not reliably give a real
# calendar date either. Convert to Date immediately (as is already done
# for nintendo_cal's own anchor above) so everything downstream operates
# on Dates consistently.
yearly_nintendo_dates <- as.Date(zoo::index(nintendo_sales)[4 * (1:19)])
first_month <- format(as.Date(zoo::index(nintendo_sales)[4]), "%m") # Note that the first element is a Q4
stopifnot(all(format(yearly_nintendo_dates, "%m") == first_month))
yearly_nintendo_mat <- zoo::coredata(nintendo_sales)[4 * (1:19), c("wii", "3ds")]

# Build an annual idx_series/idx_calendar anchored on the true first
# selected date (a Q3 endpoint), not a fabricated 1 January date.
yearly_nintendo_idx <- idx_series(yearly_nintendo_mat, start = 1L)
yearly_nintendo_cal <- idx_calendar(
  anchor = yearly_nintendo_dates[1],
  anchor_pos = 1L,
  amount = 1, unit = "years",
  posixct = TRUE
)
threeds_idx <- idx_series(idx_values(yearly_nintendo_idx)[, "3ds"], start = yearly_nintendo_idx$start)

est.start.y <- idx_to_pos(yearly_nintendo_cal, format(yearly_nintendo_dates[which(format(yearly_nintendo_dates, "%Y") == "2011")], "%Y-%m-%d"))
est.end.y   <- idx_to_pos(yearly_nintendo_cal, format(yearly_nintendo_dates[which(format(yearly_nintendo_dates, "%Y") == "2018")], "%Y-%m-%d"))

# Fit an annual dynamic Gompertz model to 3DS sales.
mod_3ds <- tsgc::SSModelDynamicGompertz(  
  Y = threeds_idx, sea.period = 0, 
  start = est.start.y, end = est.end.y, calendar = yearly_nintendo_cal
)
# Estimate the annual 3DS model and inspect the summary.
res_3ds <- tsgc::estimate(mod_3ds)
summary(res_3ds)
tsgc::print_model_diagnostics(res_3ds)

# Forecast annual 3DS log growth.
p <- tsgc::plot_log_forecast(
  res_3ds, Y = threeds_idx, n.ahead = n.forecasts,
  title = "Forecast of log growth rate of annual 3DS sales"
)
print(p)
save_plot(p, "3ds_sales_gomp_loggr_fcst.png")

# Forecast annual 3DS sales.
p <- tsgc::plot_forecast(
  res_3ds, n.ahead = n.forecasts, 
  title = "Forecast of new annual 3DS sales",
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "3ds_sales_gomp_fcst.png")

# Evaluate holdout accuracy for the annual 3DS forecast.
p <- tsgc::plot_holdout(
  res_3ds, Y = threeds_idx, n.ahead = n.forecasts,
  title = "Accuracy: Forecast of new annual 3DS sales",
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "3ds_sales_gomp_holdout.png")

#' 
#' ## 8.6 Leading Indicator Model with Annual Data: Wii to 3DS (Lead)
#' 
#' In this example, Wii is treated as the leading indicator 
#' and 3DS as the target, to illustrate the leading-indicator 
#' state-space model at annual frequency.
#' 
#' Annual lead-lag model using Wii as lead for 3DS; forecast and evaluate.
#' 

## ---- 8.6 Annual-leading: Wii to 3DS (Lead) ----
# Fit an annual leading-indicator model using Wii sales to forecast 3DS sales.
# Compute the annual Wii-to-3DS lag.
n.lag.y <- idx_to_pos(yearly_nintendo_cal, format(yearly_nintendo_dates[which(format(yearly_nintendo_dates, "%Y") == "2011")], "%Y-%m-%d")) -
  idx_to_pos(yearly_nintendo_cal, format(yearly_nintendo_dates[which(format(yearly_nintendo_dates, "%Y") == "2007")], "%Y-%m-%d"))
# Fit the annual Wii-to-3DS leading-indicator model.
mod_lead_y <- tsgc::SSModelLeadingIndicator(
  Y = yearly_nintendo_idx, sea.period = 0, n.lag = n.lag.y,
  start = est.start.y, end = est.end.y, 
  LeadIndCol = 1, calendar = yearly_nintendo_cal
)
# Estimate the annual leading-indicator model and inspect the summary.
res_lead_y <- tsgc::estimate(mod_lead_y)
summary(res_lead_y)
tsgc::print_model_diagnostics(res_lead_y)

# Forecast annual 3DS log growth from the leading-indicator model.
p <- tsgc::plot_log_forecast(
  res_lead_y, Y = yearly_nintendo_idx, n.ahead = n.forecasts, 
  title = "Forecast of log growth rate of annual 3DS sales"
)
print(p)
save_plot(p, "3ds_sales_lead_loggr_fcst.png")

# Forecast annual 3DS sales from the leading-indicator model.
p <- tsgc::plot_forecast(
  res_lead_y, n.ahead = n.forecasts, 
  title = "Forecast of new annual 3DS sales", 
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "3ds_sales_lead_fcst.png")

# Evaluate holdout accuracy for the annual leading-indicator forecast.
#
# *** SMALL-SAMPLE CAVEAT (only 2 holdout observations) ***
# n.forecasts is 2 at annual frequency here, so every accuracy/coverage
# figure below (MAPE, sMAPE, interval coverage) is computed from just 2
# holdout points. A coverage figure of 100% (or 0%) from 2 observations
# has essentially no calibration content - it says nothing reliable
# about whether the model's stated confidence level is accurate in
# general. These figures should be read as an illustration of the
# annual-frequency API, not as evidence of forecast accuracy or
# interval calibration for this series.
p <- tsgc::plot_holdout(
  res_lead_y, Y = yearly_nintendo_idx, n.ahead = n.forecasts, 
  title = paste0("Accuracy: Forecast of new annual 3DS sales\n",
                 "(illustrative only - based on n=", n.forecasts,
                 " holdout observations)"), 
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "3ds_sales_lead_holdout.png")

#' 
#' # 9. Appendix: Controlling the Plot X-Axis with `idx_axis_opts()`
#' 
#' Every plotting function used above leaves its `axis` argument unset,
#' which defaults to `mode = "auto"`: real calendar dates if the
#' calendar has `posixct = TRUE` (true of every calendar built in this
#' script), otherwise plain integer positions. `idx_axis_opts()` lets
#' you override this per plot, and applies uniformly across every
#' `plot.*`/`plot_*` function in the package. To keep the comparison
#' clear, this section reuses a single plot - the Gauteng fourteen-day
#' new-cases forecast from Section 2 - and simply varies `axis`.
#' 

## ---- 9.0 Axis showcase helper ----
# n.forecasts/plt.length are reset here to the Section 2 (Gauteng)
# values, since later sections reassign them for other frequencies.
n.forecasts <- n.forecasts.default
plt.length  <- plt.length.default

# A small wrapper so each axis mode below re-renders the same underlying
# forecast plot, varying only the `axis` argument.
showcase_axis_plot <- function(axis = NULL) {
  tsgc::plot_forecast(
    res = res_q, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
    plt.start = tail(res_q$index, 1) - plt.length,
    series.name = "cases", axis = axis
  )
}

## ---- 9.1 mode = "date" (the default here) ----
# Real calendar dates, since gauteng_cal$posixct = TRUE.
p <- showcase_axis_plot(tsgc::idx_axis_opts(mode = "date"))
save_plot(p, "axis_mode_date.png")

## ---- 9.2 mode = "position" ----
# Raw integer idx_series positions, with no calendar translation at all.
p <- showcase_axis_plot(tsgc::idx_axis_opts(mode = "position"))
save_plot(p, "axis_mode_position.png")

## ---- 9.3 mode = "steps" ----
# Steps from the calendar's anchor (here, days since the anchor position).
p <- showcase_axis_plot(tsgc::idx_axis_opts(mode = "steps"))
save_plot(p, "axis_mode_steps.png")

## ---- 9.4 mode = "time_since" ----
# The pattern-weighted calendar offset from the anchor, expressed in
# calendar$unit's rather than a raw step count - here "days since first
# recorded case". For a daily series with no gaps this is numerically
# identical to "steps", but it diverges for series with an irregular
# pattern (e.g. business days) or a non-day unit (e.g. quarters, months).
p <- showcase_axis_plot(tsgc::idx_axis_opts(mode = "time_since"))
save_plot(p, "axis_mode_time_since.png")

## ---- 9.5 Adding an info box ----
# Setting info_box = TRUE (on any of the modes above) appends a caption
# summarising the calendar: the anchor (and its name, if set), the step
# size/unit, and the pattern.
p <- showcase_axis_plot(tsgc::idx_axis_opts(mode = "steps", info_box = TRUE))
save_plot(p, "axis_mode_steps_infobox.png")

# pattern_n truncates a long pattern to its first n values in the info
# box caption; by default the full pattern is shown. Not very useful for
# Gauteng's plain daily pattern, but essential for keeping the caption
# readable with, e.g., a business-day calendar.
p <- showcase_axis_plot(
  tsgc::idx_axis_opts(mode = "steps", info_box = TRUE, pattern_n = 3)
)
save_plot(p, "axis_mode_steps_infobox_truncated.png")

#' 
#' # 10. Appendix: Compound and Heterogeneous Calendar Steps
#' 
#' Every calendar used in Sections 2-8 is built with the primary
#' `idx_calendar()` constructor, whose step size is a single
#' `amount`/`unit` pair (e.g. `amount = 1, unit = "days"`, or
#' `amount = 1, unit = "quarters"`). This covers the common case where
#' the step between consecutive `idx_series` positions is a whole
#' multiple of one calendar unit. Some series, however, are genuinely
#' spaced by a combination of units that a single `amount`/`unit` pair
#' cannot express, or by a cycle of different *kinds* of steps
#' altogether. The two subsections below illustrate these cases, plus
#' the non-calendar anchor case handled by `idx_offset_to_pos()`.
#' 

#' 
#' ## 10.1 `idx_step()` and `idx_calendar_step()`
#' 
#' `idx_step()` builds a compound step out of any combination of
#' years/quarters/months/weeks/days/hours/minutes/seconds - for example,
#' a reporting system that timestamps quarterly figures a few seconds
#' after the quarter boundary, so the step is "one quarter plus three
#' seconds". `idx_calendar_step()` is the `idx_calendar` constructor
#' that takes such a step.
#' 

## ---- 10.1 Compound step: idx_step() / idx_calendar_step() ----
qtr_step <- tsgc::idx_step(quarters = 1, seconds = 3)
qtr_step

cal_step <- tsgc::idx_calendar_step(
  anchor = as.POSIXct("2024-01-01 00:00:03", tz = "UTC"),
  step = qtr_step,
  posixct = TRUE
)
tsgc::idx_to_date(cal_step, 1:4)

# idx_to_pos() inverts this exactly as it does for the primary
# constructor. A compound step has no single amount to divide by, so
# the position is located by a binary search on the number of steps from
# the anchor, as it is for idx_calendar() calendars with numeric patterns.
tsgc::idx_to_pos(cal_step, tsgc::idx_to_date(cal_step, 4))

# Any component of an idx_step left at its default of 0 is simply
# omitted from the step, so idx_step(days = 1) is equivalent to the
# amount = 1, unit = "days" pairs used throughout this script.
# idx_step_add() is the lower-level function that applies a given
# number of whole steps to an anchor date directly; it is what an
# idx_calendar_step()-built calendar uses internally, and is not
# usually called directly by users of the package.

#' 
#' ## 10.2 `multi_step_pattern()` and `idx_calendar_multi_step()`
#' 
#' The `pattern` argument of `idx_calendar()`/`idx_calendar_step()`
#' handles a repeating cycle of *different multiples* of one step - e.g.
#' a business-day calendar, where four single-day steps (Mon-Thu) are
#' followed by a three-day step over the weekend. That still assumes
#' every step in the cycle is the same *kind* of step, just scaled
#' differently. Some data instead cycle through genuinely different
#' kinds of steps - for example, "three days, three days, then one
#' month, repeating". `multi_step_pattern()` describes this kind of
#' cycle as a sequence of `idx_step` objects, and
#' `idx_calendar_multi_step()` is the `idx_calendar` constructor that
#' uses it.
#' 

## ---- 10.2 Heterogeneous cycle: multi_step_pattern() / idx_calendar_multi_step() ----
cal_multi <- tsgc::idx_calendar_multi_step(
  anchor = as.Date("2024-01-01"),
  multi_step = tsgc::multi_step_pattern(
    tsgc::idx_step(days = 3), tsgc::idx_step(days = 3), tsgc::idx_step(months = 1)
  ),
  posixct = TRUE
)
tsgc::idx_to_date(cal_multi, 1:6)

# Because the slots in a multi_step_pattern can mix calendar-relative
# (e.g. months) and fixed-duration (e.g. days) components, there is no
# closed-form way to collapse the whole cycle into a single numeric
# offset. idx_to_date() instead walks one slot at a time from the
# anchor position, and idx_to_pos()'s inverse walks the same cycle in
# reverse.
tsgc::idx_to_pos(cal_multi, "2024-02-07")

# This makes idx_calendar_multi_step() calendars more expensive to
# convert than the other two constructors - each position takes time
# proportional to its distance from the anchor, rather than constant
# time - so it is best reserved for genuinely irregular step sequences
# that the pattern argument of the other two constructors cannot
# express. The "time_since" axis mode from Section 9 is also
# unavailable for this calendar: a multi_step_pattern has no single
# amount to express an offset in, so mode = "steps" or mode = "date"
# should be used instead.

#' 
#' ## 10.3 Non-calendar `idx_calendar` anchors and `idx_offset_to_pos()`
#' 
#' Not every `idx_series` needs a genuine calendar interpretation. An
#' `idx_calendar`'s `anchor` can be any single reference point, not just
#' a `Date`/`POSIXct` - for example, a plain numeric count of
#' picoseconds since the start of an experiment. `idx_offset_to_pos()`
#' is the inverse operation for this case: given a value already
#' expressed in the anchor's own units, it returns the integer
#' `idx_series` position that corresponds to it - the counterpart of
#' `idx_to_pos()`, which instead expects a calendar date.
#' 

## ---- 10.3 Non-calendar anchor: idx_offset_to_pos() ----
cal_ps <- tsgc::idx_calendar(
  anchor = 0, anchor_pos = 1L, amount = 2.5, unit = "picoseconds"
)
tsgc::idx_offset_to_pos(cal_ps, 12.5)

# None of the calendars used elsewhere in this script need
# idx_offset_to_pos(), since every series used above has a genuine
# Date/POSIXct calendar and so uses idx_to_pos() instead; it is
# included here for completeness, for applications - e.g. raw
# instrumentation data with no calendar meaning - where idx_series
# positions are still worth translating to and from the anchor's own
# numeric scale.

# Wrap-up
# Validate every figure written to images_dir this run: must exist,
# be non-empty, decode as a PNG, and exceed a plausible minimum size.
if (SAVE_PLOTS) {
  saved_figs <- list.files(images_dir, pattern = "\\.png$", full.names = TRUE)
  bad <- character(0)
  for (f in saved_figs) {
    ok <- tryCatch({ validate_saved_figure(f); TRUE }, error = function(e) {
      message("INVALID FIGURE: ", f, " -- ", conditionMessage(e)); FALSE
    })
    if (!ok) bad <- c(bad, f)
  }
  if (length(bad) > 0) {
    stop(length(bad), " figure(s) failed validation: ", paste(basename(bad), collapse = ", "))
  }
  message("All ", length(saved_figs), " saved figures validated (exist, non-empty, decodable).")
}
message("=== Run completed. Check 'results/Tables' and 'results/Images'. ===")