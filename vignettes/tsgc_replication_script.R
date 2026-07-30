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
#    - Supplying future xpred (xts, CSV)
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
# 7. Comparing Leading Indicator and Gompertz Models (UK–Italy)
#    7.1 Case 1: First peak window (2020-02-25 to 2020-04-01)
#        - UK-only Gompertz vs Italy→UK leading indicator
#    7.2 Cross-validation (same window as Case 1)
#        - Gompertz vs leading-indicator models with lags
#    7.3 Case 2: Extended window (2020-02-25 to 2020-04-15)
#        - Re-comparison under longer sample
#
# 8. Extensions to Other Frequencies
#    8.1 Quarterly: Wii sales (Gompertz)
#    8.2 Quarterly: Wii→Switch (leading indicator)
#    8.3 Monthly: Plus500 (Gompertz)
#    8.4 Monthly: DEGIRO→AvaTrade (leading indicator)
#    8.5 Annual: 3DS (Gompertz)
#    8.6 Annual: Wii→3DS (leading indicator)
#
# 9. Limitations, Diagnostics & Exported-File Notes
#    - Illustrative vs retrospective vs operational forecasts
#    - Cross-validation scope and interval calibration
#    - Convergence/residual diagnostics checklist
#    - Exported CSV file labelling
#
# ==========================================

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
SAVE_PLOTS   <- TRUE
SAVE_TABLES  <- TRUE
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
gen_int <- 4   # Assumed generation interval (days between infection in one case and
# secondary infections). R_t = exp(g_y,t * gen_int), where g_y,t is the
# estimated daily log-growth rate. This is an assumption, not an estimate;
# results should be checked for sensitivity to this value (see below).
ndays   <- 7   # Number of most recent daily R_t estimates to return/plot.
# NB. estimate_r0() computes R_t for every day in the sample and then
# takes tail(..., ndays): this truncates to the last `ndays` days, it
# does NOT smooth or average the growth signal. Any noise reduction in
# the daily case data happens upstream, in the Gompertz model's own
# smoothing of the growth-rate state -- not from this parameter.


# Default estimation window used in Section 2
# This is the baseline Gauteng estimation window; later holdout windows are derived from these dates.
est.start.1 <- as.Date("2021-02-01")
est.end.1   <- as.Date("2021-05-03")

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
  "tsgc","KFAS","dplyr","ggplot2","ggthemes","ggfortify","ggforce",
  "magrittr","zoo","latex2exp","xts","gridExtra","here","timetk",
  "tidyr","abind","scales","grid"
)
invisible(lapply(libs, safe_library))

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

# ---- 1.5 Plot/Save Helpers ----
# Define wrappers for saving plots. These centralise image dimensions, resolution, and the SAVE_PLOTS toggle.
# Save a ggplot with common width, height, and DPI settings when SAVE_PLOTS is TRUE.
safe_ggsave <- function(plot, filename, width = FIG_WIDTH, 
                        height = FIG_HEIGHT, dpi = FIG_DPI) {
  if (!SAVE_PLOTS) return(invisible(TRUE))
  ggplot2::ggsave(filename = filename, plot = plot, 
                  width = width, height = height, dpi = dpi)
  message("Saved plot: ", normalizePath(filename, winslash = "/"))
}

# Print each plot in interactive/rendered output and optionally save it to the Images folder.
save_plot <- function(p, fname = NULL) {
  print(p)
  if (!is.null(fname) && inherits(p, "ggplot")) 
    safe_ggsave(p, file.path(images_dir, fname))
  invisible(p)
}

# ---- 1.6 CSV Export Helpers ----
# Define CSV-export helpers for forecasts, filtered states, growth rates, R_t outputs, and the manifest used to document exported files.
# Convert time indexes from Date, yearmon, or yearqtr forms into CSV-friendly labels.
format_csv_date <- function(index_vec) {
  if (inherits(index_vec, "Date")) {
    format(index_vec, "%Y-%m-%d")
  } else {
    as.character(index_vec)
  }
}

# Turn the confidence level into a suffix, for example 68pct, for unambiguous exported column names.
confidence_suffix <- function(confidence.level) {
  as.character(round(confidence.level * 100))
}

# Export an xts object with its time index restored as an explicit Date column.
write_xts_csv <- function(x, file, columns) {
  out <- data.frame(
    Date = format_csv_date(zoo::index(x)),
    zoo::coredata(x),
    check.names = FALSE
  )
  names(out) <- c("Date", columns)
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
  
  y.hat.diff <- res$predict_level(
    n.ahead = n.ahead,
    confidence.level = confidence.level,
    sea.on = TRUE
  )
  write_xts_csv(
    y.hat.diff[, 1:3],
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
  delta.std.err <- sqrt(P.t.t[idx.level, idx.level, ])
  gamma.std.err <- sqrt(P.t.t[idx.slope, idx.slope, ])
  
  delta <- xts::xts(
    cbind(
      delta_log_growth_level = as.numeric(filtered.level),
      delta_std_error = as.numeric(delta.std.err)
    ),
    order.by = zoo::index(filtered.level)
  )
  write_xts_csv(
    delta,
    file.path(res.dir, paste0(model_slug, "_delta_filtered.csv")),
    c("delta_log_growth_level", "delta_std_error")
  )
  
  gamma <- xts::xts(
    cbind(
      gamma_trend_slope = as.numeric(filtered.slope),
      gamma_std_error = as.numeric(gamma.std.err)
    ),
    order.by = zoo::index(filtered.slope)
  )
  write_xts_csv(
    gamma,
    file.path(res.dir, paste0(model_slug, "_gamma_filtered.csv")),
    c("gamma_trend_slope", "gamma_std_error")
  )
  
  fitted.growth <- exp(filtered.level) + filtered.slope
  ci.offset <- stats::qnorm((1 - confidence.level) / 2) *
    as.numeric(gamma.std.err) %o% c(1, -1)
  growth.ci <- xts::xts(
    cbind(
      fitted_incidence_growth_rate = as.numeric(fitted.growth),
      lower = as.numeric(as.numeric(fitted.growth) + ci.offset[, 1]),
      upper = as.numeric(as.numeric(fitted.growth) + ci.offset[, 2])
    ),
    order.by = zoo::index(filtered.level)
  )
  write_xts_csv(
    growth.ci,
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

# Convenience for date windows used in plots
# Helper for plot windows: move back k periods from the last available date.
tail_date_minus <- function(index_vec, k)
  if (length(index_vec)) tail(index_vec, 1) - k else NA

# Takes the last date and goes back k days; 
# if the vector is empty, returns NA.

# ---- 1.7 knitr Defaults ----
# Configure knitr chunk defaults when the script is rendered as a vignette or report.
knitr::opts_chunk$set(
  echo       = TRUE,
  message    = TRUE,
  warning    = FALSE,
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
#' Load the Gauteng dataset and extract the cumulative cases series for modelling.
#' 

## ---- 2.1 Data -----------------
# Load the Gauteng COVID-19 example and select cumulative cases, which are the response series for the baseline Gompertz model.
# Load the built-in Gauteng dataset from tsgc.
data(gauteng, package = "tsgc")
# Select cumulative cases; the Gompertz model is fitted to cumulative counts.
cumulative_cases <- gauteng[, 1]

#' 
#' ## 2.2 Quick Inspection of Data
#' 
#' Visual sanity check of the series and index.
#' 

## ---- 2.2 Data Inspection ------------------------------------------
# Build a quick model object and plot the raw series as a visual check before formal estimation.
# Construct a model object for plotting the series before formal estimation.
mod1 <- tsgc::SSModelDynamicGompertz(Y = cumulative_cases)
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
  Y = cumulative_cases, start.date = est.start.1, 
  end.date = est.end.1
)
# Estimate the free-q model and print its summary for parameter/state diagnostics.
res_free <- tsgc::estimate(model_free); summary(res_free)

# 2.3b Diffuse prior with AR(1)
# Specify an AR(1) slope variant to impose smoother slope evolution.
model_ar1 <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, ar1 = TRUE, start.date = est.start.1, 
  end.date = est.end.1
)
# Estimate the AR(1) variant and inspect the summary.
res_ar1 <- tsgc::estimate(model_ar1); summary(res_ar1)

# 2.3c Fixed q
# Specify the preferred fixed-q model using the shared q.default value.
model_q <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, q = q.default, start.date = est.start.1, 
  end.date = est.end.1
)
# Estimate the fixed-q model used for the baseline forecasts.
res_q <- tsgc::estimate(model_q); summary(res_q)

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
  plt.start.date = tail_date_minus(res_q$index, plt.length),
  title = "Forecast of log growth rate of cases\n14-days (Gauteng)"
); print(p)

# 2.4b New cases forecast
# Forecast daily new cases in the original data scale with prediction intervals.
p <- tsgc::plot_forecast(
  res_q, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_q$index, plt.length),
  title = "Forecast of new cases\n14-days (Gauteng)", 
  series.name = "Cases"
); print(p)

# 2.4c Holdout accuracy: two weeks prior to end of sample
# Define a holdout estimation end 14 days before est.end.1
# Define a truncated estimation end date so the final 14 days can be held out for validation.
est.end.holdout <- est.end.1 - n.forecasts

# Refit ONLY for holdout evaluation on the truncated window
# Refit the fixed-q Gompertz model using only data available before the holdout period.
model_q_holdout <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases,  
  q = q.default, 
  start.date = est.start.1,
  end.date   = est.end.holdout
)
# Estimate the holdout model before comparing its forecasts with the withheld observations.
res_q_holdout <- tsgc::estimate(model_q_holdout); 
summary(res_q_holdout)

# 2.4c Holdout accuracy plot
# Plot holdout accuracy by comparing forecasts with observed values after the truncated end date.
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

# Load built-in daily Gauteng weather data for use as exogenous regressors.
data(gauteng_weather_2021, package = "tsgc")

# Subset to the estimation window [est.start.1, est.end.1]
# Restrict weather regressors to the same estimation window as the response series.
gauteng_weather_est <- get_timeframe(
  gauteng_weather_2021[, c(1, 3)],
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
  start.date = est.start.1,
  end.date   = est.end.1
)

# Estimate the weather-augmented model and inspect its fitted output.
res_weather <- tsgc::estimate(model_weather)
summary(res_weather)

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

# Example: take future rows from `gauteng_weather_2021` for the
# forecast horizon of length n.forecasts:
# Extract the future weather path needed for out-of-sample forecasts.
gauteng_weather_future <- get_timeframe(
  gauteng_weather_2021[, c(1, 3)],
  est.end.1 + 1,
  est.end.1 + n.forecasts
)

# Supply these future regressors to the fitted model:
# Attach future regressor values to the fitted model before forecasting.
tsgc::supply_xpred.new(res_weather, gauteng_weather_future)

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

# Convert the data frame into an xts object: the Date column is the index,
# and the remaining columns are the regressors supplied to the model.
# Convert the CSV data frame into xts format expected by tsgc.
gauteng_weather_future_xts <- xts(
  gauteng_weather_future_csv[, -1],
  order.by = gauteng_weather_future_csv$Date
)

# Supply CSV-based future xpred data to the model
# Supply the CSV-derived future regressors to the fitted model.
tsgc::supply_xpred.new(res_weather, gauteng_weather_future_xts)

# Once xpred has been supplied, the fitted model will generate forecasts
# that are conditional on these external regressors.

# ----------------------------------------------------
# 3.1.4 Forecasts and Accuracy Plots (with Regressors)
# ----------------------------------------------------

# Forecast log growth from the weather-augmented model.
p <- tsgc::plot_log_forecast(
  res_weather,
  Y              = cumulative_cases,
  n.ahead        = n.forecasts,
  plt.start.date = tail_date_minus(res_weather$index, plt.length),
  title          = "Forecast of log growth rate of cases\n(with regressors: weather)"
)
print(p)

# Forecast new cases conditional on the supplied future weather values.
p <- tsgc::plot_forecast(
  res_weather,
  n.ahead          = n.forecasts,
  confidence.level = CONF_LEVEL,
  plt.start.date   = tail_date_minus(res_weather$index, plt.length),
  title            = "Forecast of new cases\nwith regressors (weather), Gauteng",
  series.name      = "Cases"
)
print(p)

# Holdout accuracy: two weeks prior to end of sample
# Define a holdout estimation end 14 days before est.end.1
# Reuse the same 14-day holdout design for the model with regressors.
est.end.holdout <- est.end.1 - n.forecasts

# Refit ONLY for holdout evaluation on the truncated window
# Refit the weather-augmented model on the truncated estimation sample.
model_q_xpred_holdout <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases,  
  xpred = gauteng_weather_est,
  q = q.default, 
  start.date = est.start.1,
  end.date   = est.end.holdout
)
# Estimate the truncated model before supplying holdout-period weather values.
res_qxpred_holdout <- tsgc::estimate(model_q_xpred_holdout)

# Extract weather values covering exactly the validation horizon.
gauteng_weather_holdout <- get_timeframe(
  gauteng_weather_2021[, c(1, 3)],
  est.end.holdout + 1,
  est.end.holdout + n.forecasts
)

# Supply the holdout-period regressors required for conditional validation forecasts.
tsgc::supply_xpred.new(res_qxpred_holdout, gauteng_weather_holdout)

# Plot holdout accuracy for the model with weather regressors.
p <- tsgc::plot_holdout(
  res_qxpred_holdout,
  Y                = cumulative_cases,
  n.ahead          = n.forecasts,
  confidence.level = CONF_LEVEL,
  title            = "Accuracy: Forecast of new cases\nwith regressors",
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
r.t <- tsgc::estimate_r0(res_q, gen_int, ndays)

# ---- 4.2 Sensitivity of R_t to the assumed generation interval ----
# gen_int is an epidemiological assumption, not something estimated from
# these data, so R_t should be checked for sensitivity to it rather than
# reported only at a single assumed value.
gen_int_sensitivity <- c(3, 4, 5, 6)
r.t_sensitivity <- lapply(gen_int_sensitivity, function(g) {
  out <- tsgc::estimate_r0(res_q, g, ndays)
  names(out) <- c("Date", "Rt", "Rt_lower", "Rt_upper")
  out$gen_int <- g
  out
})
r.t_sensitivity_df <- do.call(rbind, r.t_sensitivity)
print(r.t_sensitivity_df)

if (SAVE_TABLES) {
  write.csv(
    r.t_sensitivity_df, row.names = FALSE,
    file = file.path(tables_dir, "gauteng_gompertz_q005_rt_gen_int_sensitivity.csv")
  )
  message("Saved gauteng_gompertz_q005_rt_gen_int_sensitivity.csv")
}
if (SAVE_TABLES) {
  names(r.t) <- c(
    "Date", "Rt",
    paste0("Rt_lower_", confidence_suffix(CONF_LEVEL)),
    paste0("Rt_upper_", confidence_suffix(CONF_LEVEL))
  )
  # Save the R_t table for later use in reports or supplementary material.
  write.csv(r.t, row.names = FALSE, 
            file = file.path(tables_dir, "gauteng_gompertz_q005_rt.csv"))
  message("Saved gauteng_gompertz_q005_rt.csv")
}

# Recompute R_t with plotting enabled, using the same model and epidemiological assumptions.
p <- tsgc::estimate_r0(
  res_q, gen_int, ndays, show_plot = TRUE, 
  title = "Gauteng Reproduction numbers"
)
print(p)
save_plot(p, "gauteng_rt_gomp_q005_plot.png")

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
#' If so, set the reinitialisation date at the latest sign-change prior 
#' to the first \(2\sigma\) crossing, marking the potential start 
#' of a new growth phase. 
#' 

## ---- 5.1 Reinitialisation Trigger Setup ----
# Estimate a longer Gauteng model and derive diagnostics used to identify possible reinitialisation dates.
# Extend the estimation window so the model can detect later-wave dynamics.
est.end.2 <- as.Date("2021-06-25")

# Fit a baseline model over the longer window without reinitialisation.
model_rei_base <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1, end.date = est.end.2
)
# Estimate the longer-window model and inspect the summary before deriving diagnostics.
res_rei_base <- tsgc::estimate(model_rei_base)
summary(res_rei_base)

# KFS pieces from the fitted results object
# Pull out Kalman smoother outputs, state estimates, variances, and time index for diagnostic construction.
kfs      <- res_rei_base$output
alphaHat <- kfs$alphahat
Ptt_arr  <- kfs$Ptt
idx      <- res_rei_base$index

# Smoothed slope
# Extract the smoothed slope state, which is central to identifying changes in epidemic dynamics.
smthd.slpe.full <- xts::xts(alphaHat[, "slope"], order.by = idx)

# Variance of slope from Ptt 
# Extract slope variances to form uncertainty bands around the smoothed slope.
v_slope <- drop(Ptt_arr[2, 2, ])      
idx_use <- if (length(v_slope) == length(idx) - 1) idx[-1] else idx
smoothed.P.slope <- xts::xts(as.numeric(v_slope)[seq_along(idx_use)], 
                             order.by = idx_use)

# Combine slope estimates, uncertainty, and observed series into one diagnostic object.
d2 <- cbind(smthd.slpe.full, sqrt(smoothed.P.slope),
            1.5 * sqrt(smoothed.P.slope), 2 * sqrt(smoothed.P.slope))
d2.df <- data.frame(date = index(d2), coredata(d2))
colnames(d2.df) <- c("Date", "smthd.slpe", "plus.sd.smthd.slpe",
                     "plus.sd.smthd.slpe.1.5", "plus.sd.smthd.slpe.2")
d2.df <- dplyr::filter(d2.df, Date >= as.Date("2020-10-06"))

# Identify candidate dates where the slope signal and uncertainty satisfy the trigger rule.
trigger.df <- d2.df %>%
  dplyr::mutate(prev_smthd.slpe = dplyr::lag(smthd.slpe)) %>%
  dplyr::filter(smthd.slpe > plus.sd.smthd.slpe.2 & 
                  prev_smthd.slpe < plus.sd.smthd.slpe.2)

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
                                  color = "1 SE band"), 
                     linewidth = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe.1.5, 
                                  color = "1.5 SE band"), 
                     linewidth = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe.2, 
                                  color = "2 SE band"), 
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
    subtitle = "Smoothed slope with 1, 1.5, and 2 standard-error bands",
    x = "Date",
    y = "Smoothed slope",
    caption = "Thick vertical line: reset date. Dashed vertical line: 2-SE trigger date."
  ) +
  ggplot2::scale_x_date(date_breaks = "10 days") +
  ggplot2::scale_color_manual(
    name   = "Series", 
    values = c(
      "Smoothed slope" = "red",
      "1 SE band"      = "blue",
      "1.5 SE band"    = "green",
      "2 SE band"      = "black"
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
#' Re-estimate with a chosen reinit date; produce forecasts 
#' and compare to the no-reinit baseline.
#' 

## ---- 5.2 Reinitialisation Estimation & Forecasts ----
# Refit the Gauteng model with a reinitialisation date and compare forecasts with and without reinitialisation.
# Derive reinit.date from reinit_zero.df (Section 5.1) rather than hard-coding
# it, so the chosen restart date is traceable to the trigger diagnostic above.
if (nrow(reinit_zero.df) != 1) {
  stop("reinit_zero.df did not yield exactly one candidate reset date; inspect trigger.df/reinit_zero.df from Section 5.1.")
}
reinit.date <- reinit_zero.df$Date[1]
message("Reinitialisation date derived from trigger diagnostic: ", reinit.date)

# Fit the dynamic Gompertz model with reinitialisation activated at the selected date.
model_reinit <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1, end.date = est.end.2,
  reinit.date = reinit.date
)
# Estimate the reinitialised model and inspect its summary.
res_reinit <- tsgc::estimate(model_reinit)
summary(res_reinit)

# Forecasts after reinitialisation
# Forecast log growth after allowing the model to restart at the reinitialisation date.
p <- tsgc::plot_log_forecast(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_reinit$index, plt.length),
  title = "Forecast of log growth rate of cases\nafter reinitialisation"
); print(p)

# Forecast new cases from the reinitialised model.
p <- tsgc::plot_forecast(
  res_reinit, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_reinit$index, plt.length),
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
# Select the two England series used by the leading-indicator model.
eng <- tsgc::england[, 1:2]

# Quick plot
# Create a quick leading-indicator model to visualise the lead/target relationship with n.lag = 4.
mod2 <- tsgc::SSModelLeadingIndicator(eng, n.lag = 4)
p <- plot(
  mod2, title = "Daily COVID cases and Hospitalisations\n(England)",
  series.name.lead = "Cases", series.name.target = "Hospitalisations", 
  take.log = TRUE
)
print(p)
save_plot(p, "eng_hosp_lead_cases.png")

# Estimation window
# Define the England estimation window, plotting length, lag, and forecast horizon.
est.start.eng <- as.Date("2021-04-30")
est.end.eng   <- as.Date("2021-07-24")
plt.len.eng   <- 14
n.lag         <- 4
n.forecasts   <- 7

# Define and estimate
# Specify the England leading-indicator model for hospital admissions.
out_eng <- tsgc::SSModelLeadingIndicator(
  Y = eng, n.lag = n.lag, q = NULL, LeadIndCol = 1, sea.period = 7,
  start.date = est.start.eng, end.date = est.end.eng
)
# Estimate the England leading-indicator model and inspect the fitted summary.
res_eng <- tsgc::estimate(out_eng)
summary(res_eng)

# Forecasts
# Forecast the target-series log growth rate.
p <- tsgc::plot_log_forecast(
  res_eng, Y = eng, n.ahead = n.forecasts, 
  plt.start.date = est.end.eng - plt.len.eng,
  title = "Forecast of log growth rate of hospital admissions\n(England)"
); print(p)

# Forecast hospital admissions in the original scale.
p <- tsgc::plot_forecast(
  res_eng, n.ahead = n.forecasts, 
  plt.start.date = est.end.eng - plt.len.eng,
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
  write_results_clear(
    # Export the England leading-indicator forecast and filtered-state outputs as CSV files.
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
# Prepare weather regressors for both the lead and target equations.
xpred_lead <- xpred_targ <- england_weather_2021[, 1:4]
# Fit the England leading-indicator model with xpred regressors.
mod_eng_x <- tsgc::SSModelLeadingIndicator(
  eng, n.lag = 4, xpred_lead = xpred_lead, xpred_targ = xpred_targ,
  start.date = est.start.eng, end.date = est.end.eng
)
# Estimate the weather-augmented England model and inspect the summary.
res_eng_x <- tsgc::estimate(mod_eng_x)
summary(res_eng_x)

# Supply future regressors for the lead equation.
tsgc::supply_xpred.new(res_eng_x, 
                       england_weather_2021[, 1:4], idx = "lead")
# Supply future regressors for the target equation.
tsgc::supply_xpred.new(res_eng_x, 
                       england_weather_2021[, 1:4], idx = "targ")

# Forecast log growth from the weather-augmented leading-indicator model.
p <- tsgc::plot_log_forecast(
  res_eng_x, Y = eng, n.ahead = n.forecasts, 
  plt.start.date = est.end.eng - plt.len.eng,
  title = "Forecast of log growth rate of hospital admissions\nwith regressors, England"
); print(p)

# Forecast hospital admissions with weather regressors included.
p <- tsgc::plot_forecast(
  res_eng_x, n.ahead = n.forecasts, 
  plt.start.date = est.end.eng - plt.len.eng,
  title = "Forecast of hospital admissions\nwith regressors, England", 
  series.name = "Hospital admissions"
); print(p)

# Evaluate holdout accuracy for the weather-augmented England model.
p <- tsgc::plot_holdout(
  res_eng_x, Y = eng, n.ahead = n.forecasts,
  title = "Accuracy: Forecast of hospital admissions\nwith regressors, England", 
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
#' ## 7.1 Case 1: First Peak (UK vis-à-vis Italy)
#' 

# ---- 7.1 UK vis-a-vis Italy ----
# Compare a UK-only Gompertz forecast with an Italy-to-UK leading-indicator forecast during the first peak window.

# By default, column 1 is treated as the leading indicator and 
# column 2 as the target.

# Create a quick UK-Italy leading-indicator object to inspect the lead/target timing.
ukit <- tsgc::SSModelLeadingIndicator(tsgc::ukitaly, n.lag = 4)
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
est.start <- as.Date("2020-02-25")
est.end   <- as.Date("2020-04-01")
Yuk       <- tsgc::ukitaly[, "UK"]

# Estimate a UK-only dynamic Gompertz benchmark.
res_uk_gomp1 <- tsgc::estimate(
  tsgc::SSModelDynamicGompertz(
    Y = Yuk, q = q.default,
    start.date = est.start, 
    end.date   = est.end
  )
)

# Forecast UK daily cases from the Gompertz benchmark.
p <- tsgc::plot_forecast(
  res_uk_gomp1, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz)",
  plt.start.date = tail_date_minus(res_uk_gomp1$index, plt.length),
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
n.lag <- 14
# Estimate the Italy-to-UK leading-indicator model for the same window.
res_uk_lead1 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(
    Y = tsgc::ukitaly, 
    n.lag = n.lag, sea.period = 7,
    start.date = est.start, 
    end.date   = est.end
  )
)

# Forecast UK daily cases from the leading-indicator model.
p <- tsgc::plot_forecast(
  res_uk_lead1, n.ahead = n.forecasts,
  title = "Leading indicator forecast\ndaily COVID cases in UK",
  plt.start.date = est.end - 30, series.name = "UK cases"
); print(p)

# Evaluate holdout accuracy for the leading-indicator forecast.
p <- tsgc::plot_holdout(
  res_uk_lead1, Y = tsgc::ukitaly, n.ahead = n.forecasts,
  title      = "Accuracy: Leading indicator forecast\ndaily COVID cases in UK", 
  series.name = "UK cases"
); print(p)

# Plot the Gompertz and leading-indicator forecasts together for visual comparison.
tsgc::plot_compare_forecast(
  list(res_uk_gomp1, res_uk_lead1), 
  actual = tsgc::ukitaly[, "UK"]
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
# Initialise the model list that will be passed to cross-validation.
cv_models <- list()

# Model 1: Vanilla Gompertz
# Add the fixed-q Gompertz candidate to the cross-validation set.
cv_models[["Vanilla_q"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, q = q.default, start.date = est.start, 
  end.date = est.end
)

# Model 2: Vanilla Gompertz with AR(1)
# Add the AR(1) Gompertz candidate to the cross-validation set.
cv_models[["Vanilla_ar1"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, start.date = est.start, 
  end.date = est.end, ar1 = TRUE
)

# Model 3–6: Leading Indicator with different n.lags, from 1-21
for (i in 1:21) {
  # Add leading-indicator candidates over a grid of lag choices.
  cv_models[[paste0("Lag", i)]] <- tsgc::SSModelLeadingIndicator(
    Y = tsgc::ukitaly, start.date = est.start, 
    end.date = est.end, n.lag = i
  )
}

# Run cross-validation (Case 1 window)
# Run rolling/holdout cross-validation and score models using sMAPE over a 14-step horizon.
tsgc::cross_val(
  Y           = tsgc::ukitaly, 
  model_list  = cv_models, 
  est.end.date = est.end,
  criterion    = "smape",
  n.ahead      = 14,
  n.estimate   = 5, 
  gap          = 2
)

# -------------------------------
# 7.3 Case 2: Extended Window
# -------------------------------

#' 
#' ## 7.3 Case 2: Extended Window (UK vis-à-vis Italy)
#' 
#' Re-estimate both models on an extended sample window and 
#' compare forecasts and accuracy.
#' 

# Case 2: Extended window
# Case 2 extends the estimation window to 15 April 2020 to test sensitivity to a longer sample.
est.start <- as.Date("2020-02-25")
est.end   <- as.Date("2020-04-15")

# Estimate the extended-window UK Gompertz model.
res_uk_gomp2 <- tsgc::estimate(
  tsgc::SSModelDynamicGompertz(
    Y = Yuk, q = q.default,
    start.date = est.start, end.date = est.end
  )
)

# Forecast UK cases from the extended-window Gompertz model.
p <- tsgc::plot_forecast(
  res_uk_gomp2, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz, extended)",
  plt.start.date = tail_date_minus(res_uk_gomp2$index, plt.length),
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
res_uk_lead2 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(
    Y = tsgc::ukitaly, n.lag = 14,
    start.date = est.start, end.date = est.end
  )
)

# Forecast UK cases from the extended-window leading-indicator model.
p <- tsgc::plot_forecast(
  res_uk_lead2, n.ahead = n.forecasts,
  title = "Forecast of daily COVID cases\nUK (Leading indicator model, extended)",
  plt.start.date = est.end - plt.length, 
  series.name = "UK cases"
); print(p)

# Evaluate holdout accuracy for the extended-window leading-indicator model.
p <- tsgc::plot_holdout(
  res_uk_lead2, Y = tsgc::ukitaly, n.ahead = n.forecasts,
  title      = "Accuracy: Forecast of daily COVID cases\nUK (Leading indicator model, extended)", 
  series.name = "UK cases"
); print(p)

# Compare the two extended-window forecasts directly.
tsgc::plot_compare_forecast(
  list(res_uk_gomp2, res_uk_lead2), 
  actual = tsgc::ukitaly[, "UK"]
)

# ==========================
# 8. Other Data Frequencies
# ==========================

#' 
#' # 8. Extensions to Other Data Frequencies 
#' 
#' Demonstrate quarterly, monthly, and annual use-cases 
#' with appropriate seasonal settings and lead–lag structures.
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
# Load Nintendo quarterly sales data and select Wii sales for the Gompertz example.
data(nintendo_sales, package = "tsgc")
wii <- nintendo_sales[, 1]

# Note: The tsgc function requires the input cumulative series to be strictly 
# increasing in time. If the cumulative values exhibit plateaus—as in 
# the case of the Wii series—it is necessary to add small increments 
# to eliminate flat segments and allow model estimation. Note that a 
# strictly increasing cumulative series also implies that the underlying 
# (non-cumulative) series must be strictly positive.

# Model estimated for a strictly increasing segment
# Use a four-quarter forecast horizon and yearqtr dates for the quarterly model.
n.forecasts <- 4
est.start.q <- zoo::as.yearqtr("2006 Q4")
est.end.q   <- zoo::as.yearqtr("2010 Q3")

# Fit a quarterly dynamic Gompertz model to Wii sales.
mod_wii <- tsgc::SSModelDynamicGompertz(  
  Y = wii, sea.period = 4, start.date = est.start.q, 
  end.date = est.end.q
)
# Estimate the quarterly Wii model and inspect the summary.
res_wii <- tsgc::estimate(mod_wii)
summary(res_wii)

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
p <- tsgc::plot_holdout(
  res_wii, Y = wii, n.ahead = n.forecasts, 
  title = "Accuracy: Forecast of new Wii sales",
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "wii_sales_gomp_holdout.png")

#' ## 8.2 Leading Indicator Model with Quarterly Data: Wii to Switch
#' 
#'  Nintendo Wii was launched in 2006, and the Nintendo Switch in 2017.
#' 
#' Model quarterly lead–lag between Wii and Switch; forecast and evaluate.
#' 

## ----  8.2 Quarterly: Wii to Switch (Lead) ----
# Use Wii sales as a quarterly leading indicator for Switch sales and forecast the target series.
# Set the quarterly leading-indicator horizon and estimation window for Switch sales.
n.forecasts   <- 8
est.start.q2  <- zoo::as.yearqtr("2017 Q1")
est.end.q2    <- zoo::as.yearqtr("2019 Q4")
# Compute the lead-lag distance between Wii and Switch launches in quarters.
n.lag.q       <- round((zoo::as.yearqtr("2017 Q1") - zoo::as.yearqtr("2006 Q4"))*4)

# Column 1 (Wii) is treated as the lead; column 2 (Switch) as the target.
# Select the lead and target quarterly sales series.
y_q <- nintendo_sales[, c("wii", "switch_all")]
# Fit the quarterly Wii-to-Switch leading-indicator model.
mod_switch <- tsgc::SSModelLeadingIndicator(
  Y = y_q, sea.period = 4, n.lag = n.lag.q, 
  start.date = est.start.q2, end.date = est.end.q2
)
# Estimate the quarterly leading-indicator model and inspect the summary.
res_switch <- tsgc::estimate(mod_switch)
summary(res_switch)

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
p <- tsgc::plot_holdout(
  res_switch, Y = y_q, n.ahead = n.forecasts, 
  title = "Accuracy: Forecast of new Switch sales",
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
# Load monthly e-trading app downloads and select Plus500 for the Gompertz example.
data(etrading_apps, package = "tsgc")
Plus500 <- etrading_apps[, 1]

# Use a four-month horizon and yearmon dates for the monthly model.
n.forecasts <- 4
est.start.m <- zoo::as.yearmon(2016)
est.end.m   <- zoo::as.yearmon(2021)

# Fit a monthly dynamic Gompertz model to Plus500 downloads.
mod_500 <- tsgc::SSModelDynamicGompertz(
  Y = Plus500, sea.period = 12, start.date = est.start.m, 
  end.date = est.end.m
)
# Estimate the monthly Plus500 model and inspect the summary.
res_500 <- tsgc::estimate(mod_500)
summary(res_500)

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
p <- tsgc::plot_holdout(
  res_500, Y = Plus500, n.ahead = n.forecasts, 
  title = "Accuracy: Forecast of new Plus500 downloads",
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
# Set the monthly leading-indicator horizon, dates, and lag for AvaTrade.
n.forecasts  <- 4
est.start.m2 <- zoo::as.yearmon(2017.5)
est.end.m2   <- zoo::as.yearmon(2021 + 1/12)
# Compute the DEGIRO-to-AvaTrade lag in months.
n.lag.m      <- round((zoo::as.yearmon(2017.5) - zoo::as.yearmon(2017))*12)

# Select the monthly lead and target app-download series.
y_m <- etrading_apps[, c("DEGIRO", "AvaTrade")]
# Fit the monthly DEGIRO-to-AvaTrade leading-indicator model.
mod_500_lead <- tsgc::SSModelLeadingIndicator(
  Y = y_m, sea.period = 12, n.lag = n.lag.m, 
  start.date = est.start.m2, end.date = est.end.m2
)
# Estimate the monthly leading-indicator model and inspect the summary.
res_500_lead <- tsgc::estimate(mod_500_lead)
summary(res_500_lead)

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
p <- tsgc::plot_holdout(
  res_500_lead, Y = y_m, n.ahead = n.forecasts, 
  title = "Accuracy: Forecast of new AvaTrade downloads", 
  series.name = "Downloads"
)
print(p)
save_plot(p, "avatrade_downloads_lead_holdout.png")

#' 
#' ## 8.5 Annual: 3DS
#' 
#' 3DS is Nintendo's handheld console (released 2011).
#' Here we convert quarterly global sales to annual frequency 
#' to demonstrate annual Gompertz modelling.
#' 
#' Annual modelling with `sea.period = 0`; 
#' forecast and evaluate 3DS annual series.
#' 

## ---- 8.5 Annual: 3DS ----
# Aggregate quarterly Nintendo data to annual frequency and fit an annual Gompertz model.
# Use a two-year forecast horizon and annualised yearmon dates for the annual model.
n.forecasts <- 2
# NB. dates below use December (Q4) to match the year-end dating of the
# annual series constructed below (yearly_nintendo_dates); as.yearmon(2011)
# would default to January and misalign with the December-dated observations.
est.start.y <- zoo::as.yearmon("2011-12")
est.end.y   <- zoo::as.yearmon("2018-12")

# Convert quarterly to yearly: take the Q4 (year-end) observation of each year.
# nintendo_sales holds CUMULATIVE sales, so a calendar-year snapshot is the
# year-end level, not a sum of quarters (summing cumulative values would
# double-count).
# Row indices are derived from the actual row labels (e.g. "2005 Q4") rather
# than assumed from a hard-coded start quarter. Confirmed: nintendo_sales
# starts at 2004 Q4, not 2005 Q1 -- so `4 * (1:19)` selects 2005 Q3, 2006 Q3,
# ..., i.e. Q3-to-Q3 values mislabeled as annual/January. Deriving the index
# from row labels avoids re-introducing that bug if the data's start quarter
# ever changes.
# Row indices are derived from the actual time index (a `yearqtr`, via
# zoo::index()) rather than assumed from a hard-coded start quarter.
# NB. nintendo_sales is an xts object, so rownames() is empty -- its labels
# live in the yearqtr index, not in character rownames. Confirmed:
# nintendo_sales starts at 2004 Q4, not 2005 Q1 -- so `4 * (1:19)` selects
# 2005 Q3, 2006 Q3, ..., i.e. Q3-to-Q3 values mislabeled as annual/January.
# Matching directly against the yearqtr index avoids re-introducing that bug
# if the data's start quarter ever changes.
target_years <- 2005:2023
nintendo_idx <- zoo::index(nintendo_sales)
target_yearqtrs <- zoo::as.yearqtr(paste0(target_years, " Q4"))
yearly_nintendo_idx <- match(target_yearqtrs, nintendo_idx)
if (any(is.na(yearly_nintendo_idx))) {
  stop("Could not find a Q4 observation for every year in 2005:2023 in nintendo_sales; ",
       "inspect zoo::index(nintendo_sales) and adjust target_years/matching.")
}
# Subsample quarterly Nintendo sales to true calendar year-end (Q4) observations.
yearly_nintendo      <- nintendo_sales[yearly_nintendo_idx, 
                                       c("wii", "3ds")]
# Label with December of each year, since these are year-end (Q4) snapshots,
# not January/annual-average values.
yearly_nintendo_dates <- zoo::as.yearmon(paste0(target_years, "-12"))
# Build annual xts series for 3DS and the wider Nintendo sales panel.
threeds_xts          <- xts::xts(
  zoo::coredata(yearly_nintendo[, "3ds"]), 
  order.by = yearly_nintendo_dates
)
yearly_nintendo_xts  <- xts::xts(
  zoo::coredata(yearly_nintendo), 
  order.by = yearly_nintendo_dates
)

# Fit an annual dynamic Gompertz model to 3DS sales.
mod_3ds <- tsgc::SSModelDynamicGompertz(  
  Y = threeds_xts, sea.period = 0, 
  start.date = est.start.y, end.date = est.end.y
)
# Estimate the annual 3DS model and inspect the summary.
res_3ds <- tsgc::estimate(mod_3ds)
summary(res_3ds)

# Forecast annual 3DS log growth.
p <- tsgc::plot_log_forecast(
  res_3ds, Y = threeds_xts, n.ahead = n.forecasts,
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
  res_3ds, Y = threeds_xts, n.ahead = n.forecasts,
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
#' Annual lead–lag model using Wii as lead for 3DS; forecast and evaluate.
#' 

## ---- 8.6 Annual-leading: Wii to 3DS (Lead) ----
# Fit an annual leading-indicator model using Wii sales to forecast 3DS sales.
# Compute the annual Wii-to-3DS lag.
n.lag.y <- zoo::as.yearmon(2011) - zoo::as.yearmon(2007)
# Fit the annual Wii-to-3DS leading-indicator model.
mod_lead_y <- tsgc::SSModelLeadingIndicator(
  Y = yearly_nintendo_xts, sea.period = 0, n.lag = n.lag.y,
  start.date = est.start.y, end.date = est.end.y, 
  LeadIndCol = 1
)
# Estimate the annual leading-indicator model and inspect the summary.
res_lead_y <- tsgc::estimate(mod_lead_y)
summary(res_lead_y)

# Forecast annual 3DS log growth from the leading-indicator model.
p <- tsgc::plot_log_forecast(
  res_lead_y, Y = yearly_nintendo_xts, n.ahead = n.forecasts, 
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
p <- tsgc::plot_holdout(
  res_lead_y, Y = yearly_nintendo_xts, n.ahead = n.forecasts, 
  title = "Accuracy: Forecast of new annual 3DS sales", 
  series.name = "Sales (million units)"
)
print(p)
save_plot(p, "3ds_sales_lead_holdout.png")

#' 
#' # 9. Limitations, Diagnostics & Exported-File Notes
#' 
#' ## 9.1 Illustrative vs. retrospective vs. operational forecasts
#' 
#' Not all examples above represent the same evidentiary standard:
#' 
#' - **Operational-style validation** (closest to a genuine real-time
#'   test): the Gauteng and England daily holdout evaluations, where the
#'   model is estimated up to a cutoff and evaluated against genuinely
#'   subsequent observations. Even here, the weather-regressor holdout
#'   uses **realised** future weather (an oracle/conditional forecast),
#'   not an archived weather forecast available at the time.
#' - **Retrospective/conditional forecasts**: any example where the
#'   "future" regressor values (e.g. weather) are the actual realised
#'   values rather than what would have been available at the forecast
#'   origin. These demonstrate the model's conditional forecasting
#'   mechanics, not its unconditional real-world accuracy.
#' - **API/illustrative demonstrations, not empirical validation**: the
#'   quarterly, monthly, and especially annual examples (Wii, Plus500,
#'   DEGIRO->AvaTrade, 3DS, Wii->3DS). Sample sizes here are small by
#'   construction -- the annual leading-indicator and Gompertz examples
#'   validate on only 2 holdout points. These sections should be read as
#'   "here is how to call the function," not "here is evidence the model
#'   works well at annual frequency."
#' 
#' ## 9.2 Cross-validation scope
#' 
#' The lag-selection cross-validation (Section 7.2) tests multiple lag
#' candidates on a small number of closely-spaced origins with overlapping
#' 14-day horizons, then selects the best lag on those same folds -- this
#' has no held-out outer evaluation and no naive/persistence benchmark. If
#' a later case study uses a different lag than the cross-validation
#' winner, that deviation should be explained explicitly; otherwise use
#' the CV-selected lag. The selection procedure should be validated
#' against an outer holdout or a nested rolling-origin scheme, with naive
#' benchmarks included, before its result is treated as a tuned
#' hyperparameter rather than a demonstration.
#' 
#' ## 9.3 Interval calibration
#' 
#' Reported coverage rates can vary widely across examples. Nominal 68%
#' intervals that are poorly calibrated in either direction (too narrow or
#' too wide, as seen when comparing across the examples in this document)
#' suggest the uncertainty quantification should not be taken at face
#' value without further diagnostic work (Section 9.4).
#' 
#' ## 9.4 Convergence, boundary, and residual diagnostics checklist
#' 
#' For each headline model above (Gauteng fixed-q, England leading
#' indicator, England + weather), before relying on its output, check:
#' 
#' - **Convergence**: did `tsgc::estimate()`'s optimiser converge without
#'   warnings? Note the number of iterations if available.
#' - **Boundary estimates**: is the estimated `q` at or near 0 (a boundary
#'   of the parameter space)? This is common in Gompertz-curve fits and
#'   affects the validity of standard-error-based inference on `q`.
#' - **Residual diagnostics**: standardised residual plots, and at
#'   minimum a check for residual autocorrelation (e.g. Ljung-Box) and
#'   heteroskedasticity.
#' - **Warnings**: warnings should be visible during model fitting and
#'   validation, and suppressed individually -- with a documented reason
#'   -- only after the underlying issue is understood, rather than
#'   suppressed globally.
#' 
#' ## 9.5 Exported file labelling
#' 
#' Files matching `*_filtered.csv` are generated via
#' `predict_all(..., return.all = TRUE)` and therefore contain both
#' filtered (in-sample) and forecast-period rows. Split these into
#' separate `*_filtered.csv` / `*_forecast.csv` outputs, or rename the
#' existing files to make the mixed content explicit (e.g.
#' `*_filtered_and_forecast.csv`), so downstream users don't assume every
#' row is an in-sample filtered estimate.
#' 
#' ## 9.6 sMAPE scale
#' 
#' Figures in this document use
#' \(\text{sMAPE} = 100 \cdot |Actual - Forecast| / (Actual + Forecast)\),
#' scaled 0-100. This is **half the scale** of the more common
#' \(200 \cdot |Actual - Forecast| / (Actual + Forecast)\) (0-200)
#' convention. This is now labelled explicitly wherever sMAPE is reported
#' (plot subtitles read "sMAPE (0-100 scale)", and `mapes()` return lists
#' include a `smape_scale` field) -- treat any sMAPE figure from this
#' package as non-comparable to externally reported sMAPE values unless
#' the scale is confirmed to match.
#' 
#' # Wrap-up
message("=== Run completed. Check 'results/Tables' and 'results/Images'. ===")