#' ---
#' title: "Replication: Time Series Growth Curves (`tsgc`)"
#' author: "Ashby, Harvey, Kattuman, Tang, Thamotheram"
#'=======================================================
#' 
#' # Setup & Utilities
#' 
#' Centralise parameters, load libraries, set global options, define paths, and provide helper functions for saving plots and controlling chunk defaults.
#' 
## ----1-0-preamble, include=TRUE----------------------------------------------------------------------------
# ===== 1. Setup & Utilities =====

# ---- 1.1 Parameters & Toggles ----
SAVE_PLOTS   <- TRUE
SAVE_TABLES  <- TRUE
FIG_WIDTH    <- 10
FIG_HEIGHT   <- 6
FIG_DPI      <- 300
CONF_LEVEL   <- 0.68

# Core analysis parameters 
n.forecasts.default <- 14
q.default           <- 0.005
plt.length.default  <- 30

# Reproduction Number (R_t) parameters
gen_int <- 4
ndays   <- 7

# Default estimation window used in Section 2
est.start.1 <- as.Date("2021-02-01")
est.end.1   <- as.Date("2021-04-19")

# ---- 1.2 Libraries (quiet require) ----
safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install before running.", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

libs <- c(
  "tsgc","KFAS","dplyr","ggplot2","ggthemes","ggfortify","ggforce",
  "magrittr","zoo","latex2exp","xts","gridExtra","here","timetk",
  "tidyr","abind","scales","grid"
)
invisible(lapply(libs, safe_library))

# ---- 1.3 Global Options & Theme ----
theme_set(ggthemes::theme_economist_white(gray_bg = FALSE, base_size = 16))
options(scipen = 7)  
set.seed(123)

# ---- 1.4 Paths & Directories ----
setwd("C:/tsgc-incorporatedClasses")
base_path   <- getwd()
results_dir <- file.path(base_path, "results")
tables_dir  <- file.path(results_dir, "Tables")
images_dir  <- file.path(results_dir, "Images")

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- tryCatch({ dir.create(path, recursive = TRUE); TRUE },
                   error = function(e) { message("Failed to create: ", path, " :: ", e$message); FALSE })
    if (!ok) stop("Could not create required directory: ", path)
  }
}
invisible(lapply(list(results_dir, tables_dir, images_dir), ensure_dir))

# ---- 1.5 Plot/Save Helpers ----
safe_ggsave <- function(plot, filename, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI) {
  if (!SAVE_PLOTS) return(invisible(TRUE))
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved plot: ", normalizePath(filename, winslash = "/"))
}
save_plot <- function(p, fname = NULL) {
  print(p)
  if (!is.null(fname) && inherits(p, "ggplot")) safe_ggsave(p, file.path(images_dir, fname))
  invisible(p)
}
# Convenience for date windows used in plots
tail_date_minus <- function(index_vec, k) if (length(index_vec)) tail(index_vec, 1) - k else NA

# ---- 1.6 knitr Defaults ----
knitr::opts_chunk$set(
  echo = TRUE, message = TRUE, warning = FALSE,
  fig.align = "center",
  fig.width = FIG_WIDTH, fig.height = FIG_HEIGHT
)

#' 
#' # Baseline Gompertz Model: Gauteng
#' 
#' The data used in this replication relates to daily COVID-19 cases in Gauteng province in
#' South Africa. We begin by fitting a baseline dynamic Gompertz model to cumulative COVID-19 cases. We then estimate variants of the baseline model, forecast and evaluate forecasts.
#' 
#' ## Data
#' 
#' Load the Gauteng dataset and extract the cumulative cases series for modelling.
#' 
## ---- 2.1 Data ----
data(gauteng, package = "tsgc")
stopifnot(ncol(gauteng) >= 1)
cumulative_cases <- gauteng[, 1]

#' 
#' ## Quick Inspection of data
#' 
#' Visual sanity check of the series and index.
#' 
## ---- 2.2 Data Inspection ----
mod1 <- tsgc::SSModelDynamicGompertz(Y = cumulative_cases)
p <- plot(mod1, title = "Gauteng daily cases", series.name = "cases")
print(p)
# save_plot(p, "gauteng_cases_MA.png")

#' ## Estimation: Gompertz model variants
#' Three Gompertz model variants: 
#' (a) diffuse/free q: flexible baseline, 
#' (b) AR(1) slope: smoother slope dynamics, 
#' (c) fixed q: chosen based on experience.
#' 
#' 
## ---- 2.3 Estimation ----
# 2.3a Diffuse prior (free q)
#model_free <- tsgc::SSModelDynamicGompertz$new(
model_free <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, start.date = est.start.1, end.date = est.end.1
)
res_free <- tsgc::estimate(model_free); summary(res_free)

# 2.3b Diffuse prior with AR(1)
#model_ar1 <- tsgc::SSModelDynamicGompertz$new(
model_ar1 <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, ar1 = TRUE, start.date = est.start.1, end.date = est.end.1
)
res_ar1 <- tsgc::estimate(model_ar1); summary(res_ar1)

# 2.3c Fixed q
# model_q <- tsgc::SSModelDynamicGompertz$new(
model_q <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases, q = q.default, start.date = est.start.1, end.date = est.end.1
)
res_q <- tsgc::estimate(model_q); summary(res_q)

#' 
#' ## Forecasts & Accuracy
#' Produce forecasts (log growth and levels) from the fixed-q model and evaluate holdout accuracy.
#' 
## ---- 2.4 Forecasts & Accuracy ----
n.forecasts <- n.forecasts.default
plt.length  <- plt.length.default

# 2.4a Log growth forecast (fixed q)
p <- tsgc::plot_log_forecast(
  res_q, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_q$index, plt.length),
  title = "Forecast of log growth rate of cases\n14-days (Gauteng)"
); print(p)

# 2.4b New cases forecast
p <- tsgc::plot_forecast(
  res_q, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_q$index, plt.length),
  title = "Forecast of new cases\n14-days (Gauteng)", series.name = "Cases"
); print(p)


# 2.4c Holdout accuracy : two weeks prior to end of sample
# Define a holdout estimation end 14 days before est.end.1
est.end.holdout <- est.end.1 - n.forecasts

# Refit ONLY for holdout evaluation on the truncated window
model_q_holdout <- tsgc::SSModelDynamicGompertz(
  Y = cumulative_cases,  
  q = q.default, 
  start.date = est.start.1,
  end.date   = est.end.holdout
)
res_q_holdout <- tsgc::estimate(model_q_holdout); summary(res_q_holdout)

# 2.4c Holdout accuracy 
p <- tsgc::plot_holdout(
  res_q_holdout, Y = cumulative_cases, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases\n14-days (Gauteng)", series.name = "PM10"
); print(p)

if (SAVE_TABLES) {
  tsgc::write_results(res = res_q, res.dir = tables_dir, 
                      n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
                      prefix= "gauteng_gomp_q005_")
  message("Saved results for: gauteng_gomp_q005")
}

#' Note: The Symmetric Mean Absolute Percentage Error (sMAPE) is a scale-free, 
#' symmetric accuracy measure that ranges from 0% to 100%. It complements MAPE, 
#' which tends to overstate forecast errors when actual values are small. 
#' However, sMAPE can also become unstable when both actual and forecast 
#' values are very small.
#' 
#' 
#' # Gauteng: Gompertz model with exogenous regressors
#' 
#' Augment Gompertz with weather regressors; re-estimate and compare forecasts and accuracy.
#' 
## ---- 3.1 Estimation with Regressors: xpred----

# Data on future values of exogenous variables (e.g., weather) 
# can be used as regressors for generating out-of-sample forecasts with tsgc.   

# In this example future xpred observations are canned in `gauteng_weather_2021`.
# We use daily weather variables as x_t (columns 1 and 3 in gauteng_weather_2021), 
# and keep the estimation window unchanged to isolate the effect of regressors.

data(gauteng_weather_2021, package = "tsgc")
gauteng_weather_est <- get_timeframe(gauteng_weather_2021[, c(1, 3)], est.start.1, est.end.1)
head(gauteng_weather_est)

# If data on future values of exogenous variables to be used as regressors 
# are to be supplied by reading in data, a csv file containing  `n.forecast` 
# observations of these regressors can be read in. The file must include 
# a 'Date' column and  one or more numeric columns. The resulting data frame 
# must be converted into an `xts` object named `xpred_new`. 
# This object input to the fitted model object will generate forecasts 
# based on the  external regressors.

model_weather <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases,
  xpred = gauteng_weather_est,
  start.date = est.start.1,
  end.date   = est.end.1
)
res_weather <- tsgc::estimate(model_weather)
summary(res_weather)

# `supply_xpred.new` is used to provide future xpred values to the model
# so that it can produce forecasts. 
#
# The object `res_weather` contains only the weather data used during the 
# estimation (training) period. To generate predictions, the model also 
# needs xpred values **beyond** this estimation window.
#
# Therefore, `supply_xpred.new(res_weather, gauteng_weather)` inserts 
# additional xpred data into the model for the forecast horizon. 
# These future xpred values can come from weather forecasts (e.g., from 
# the Met Office) or from hypothetical numbers used in best/worst-case 
# scenario analysis.
gauteng_weather_future<-get_timeframe(gauteng_weather_2021[, c(1, 3)], est.end.1+1,est.end.1+n.forecasts)
tsgc::supply_xpred.new(res_weather, gauteng_weather)

# Forecasts (xpred)
p <- tsgc::plot_log_forecast(
  res_weather, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_weather$index, plt.length),
  title = "Forecast of log growth rate of cases\n(with regressors - weather)"
); print(p)

p <- tsgc::plot_forecast(
  res_weather, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_weather$index, plt.length),
  title = "Forecast of new cases\nwith regressors (weather), Gauteng",
  series.name = "Cases"
); print(p)

p <- tsgc::plot_holdout(
  res_weather, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases\nwith regressors", series.name = "cases"
); print(p)

p <- tsgc::plot_compare_forecast(list(res_free, res_q, res_ar1, res_weather), actual = cumulative_cases)
print(p)

#' 
#' # Reproduction Number (R_t)
#' Map Gompertz estimates to effective reproduction numbers using a specified generation interval (4 days),  smoothed over a 7‑day window.
#' 
## ---- 4.1 Transform Gompertz estimates to R_t ----
r.t <- tsgc::estimate_r0(res_q, gen_int, ndays)
if (SAVE_TABLES) {
  write.csv(r.t, row.names = FALSE, file = file.path(tables_dir, "gauteng_rt_gomp_q005.csv"))
  message("Saved gauteng_rt_gomp_q005.csv")
}

p <- tsgc::estimate_r0(res_q, gen_int, ndays, show_plot = TRUE, title = "Gauteng Reproduction numbers")
print(p)
# save_plot(p, "gauteng_rt_gomp_q005_plot.png")

 
#' # Reinitialisation for subsequent waves
#' Identify a trigger for reinitialising based on smoothed slope uncertainty; 
#' re-estimate with reinitialisation and assess forecasts vs baseline.
#' 
#' ## Trigger diagnostics
#' Compute smoothed slope and uncertainty; detect threshold crossings 
#' as reinitialisation triggers. 
#' Trigger is the first point where the estimated slope rises above 
#' its own 2-sigma upper confidence band after being below it 
#' in the previous period. 
#' If so, set the re-initialisation date at the latest sign-change prior 
#' to the first \(2\sigma\) crossing. 
#' Marking the potential the start of a new growth phase. 
#' 
#' 
## ---- 5.1 Reinitialisation Trigger Setup ----
est.end.2 <- as.Date("2021-06-25")

model_rei_base <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1, end.date = est.end.2
)
res_rei_base <- tsgc::estimate(model_rei_base)
summary(res_rei_base)

# KFS pieces via package accessors
kfs      <- tsgc::output(res_rei_base)
alphaHat <- tsgc::alphahat(kfs)   
Ptt_arr  <- tsgc::Ptt(kfs)        
idx      <- res_rei_base$index

# Smoothed slope
smthd.slpe.full <- xts::xts(alphaHat[, "slope"], order.by = idx)

# Variance of slope from Ptt 
v_slope <- drop(Ptt_arr[2, 2, ])      
idx_use <- if (length(v_slope) == length(idx) - 1) idx[-1] else idx
smoothed.P.slope <- xts::xts(as.numeric(v_slope)[seq_along(idx_use)], order.by = idx_use)

d2 <- cbind(smthd.slpe.full, sqrt(smoothed.P.slope),
            1.5 * sqrt(smoothed.P.slope), 2 * sqrt(smoothed.P.slope))
d2.df <- data.frame(date = index(d2), coredata(d2))
colnames(d2.df) <- c("Date", "smthd.slpe", "plus.sd.smthd.slpe",
                     "plus.sd.smthd.slpe.1.5", "plus.sd.smthd.slpe.2")
d2.df <- dplyr::filter(d2.df, Date >= as.Date("2020-10-06"))

trigger.df <- d2.df %>%
  dplyr::mutate(prev_smthd.slpe = dplyr::lag(smthd.slpe)) %>%
  dplyr::filter(smthd.slpe > plus.sd.smthd.slpe.2 & prev_smthd.slpe < plus.sd.smthd.slpe.2)

reinit_zero.df <- d2.df %>%
  dplyr::mutate(prev_smthd.slpe = dplyr::lag(smthd.slpe)) %>%
  dplyr::filter(Date < min(trigger.df$Date) & (smthd.slpe > 0 & prev_smthd.slpe < 0)) %>%
  dplyr::arrange(dplyr::desc(Date)) %>%
  dplyr::slice(1)

#' 
## ---- 5.1 Reinitialisation Trigger Plot ----
p_trigger <-
  ggplot2::ggplot(data = d2.df[d2.df$Date > as.Date("2021-02-11"), ], ggplot2::aes(x = Date)) +
  ggplot2::geom_line(ggplot2::aes(y = smthd.slpe, color = "smthd.slpe"), linewidth = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe, color = "plus.sd.smthd.slpe"), linewidth = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe.1.5, color = "plus.sd.smthd.slpe.1.5"), linewidth = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = plus.sd.smthd.slpe.2, color = "plus.sd.smthd.slpe.2"), linewidth = 0.5) +
  ggplot2::scale_y_continuous(n.breaks = 10) +
  ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  ggplot2::geom_vline(data = trigger.df, ggplot2::aes(xintercept = Date), linewidth = 0.5, color = "black") +
  ggplot2::geom_vline(data = reinit_zero.df, ggplot2::aes(xintercept = Date), linewidth = 1, color = "black") +
  ggplot2::xlab("Day") + ggplot2::ylab("Slope") +
  ggplot2::scale_x_date(date_breaks = "10 days") +
  ggplot2::scale_color_manual(name = '', values = c('smthd.slpe' = 'red',
                                                    'plus.sd.smthd.slpe' = 'blue',
                                                    'plus.sd.smthd.slpe.1.5' = 'green',
                                                    'plus.sd.smthd.slpe.2' = 'black')) +
  ggplot2::theme_light(base_size = 11) +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 2),
                 legend.text  = ggplot2::element_text(size = 6),
                 axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                 plot.title   = ggplot2::element_text(face = "bold"))
print(p_trigger)
# save_plot(p_trigger, "gauteng_cases_gomp_q005_reinit_trigger.png")

#' 
#' ## Forecasts & Accuracy (post-reinitialisation)
#' Re-estimate with a chosen reinit date; produce forecasts 
#' and compare to the no-reinit baseline.
#' 
## ---- 5.2 Reinitialisation Estimation & Forecasts ----
# Set reinit date (could also take from trigger.df/reinit_zero.df)
reinit.date <- as.Date("2021-04-21")

model_reinit <- tsgc::SSModelDynamicGompertz(  
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1, end.date = est.end.2,
  reinit.date = reinit.date
)
res_reinit <- tsgc::estimate(model_reinit)
summary(res_reinit)

# Forecasts after reinitialistion
p <- tsgc::plot_log_forecast(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_reinit$index, plt.length),
  title = "Forecast of log growth rate of cases\nafter reinitialisation"
); print(p)

p <- tsgc::plot_forecast(
  res_reinit, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_reinit$index, plt.length),
  title = "Forecast of new cases\nafter reinitialisation", series.name = "cases"
); print(p)

p <- tsgc::plot_holdout(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL, 
  title = "Accuracy: Forecast of new cases\nwith reinitialisation", series.name = "cases"
); print(p)

# Compare holdout with/without reinitialisation (baseline = res_rei_base)
tsgc::plot_holdout(
  res_rei_base, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL, title = "Forecast without reinitialisation (baseline window)"
)
tsgc::plot_compare_forecast(list(res_rei_base, res_reinit), actual = cumulative_cases)

#' 
#' # Leading Indicator model: England (Daily)
#' Model how cases (lead) anticipate hospitalisations (target) at lag L, 
#' with weekly seasonality; produce forecasts and accuracy metrics. 
#' 
#' ## Leading Indicator: Baseline model
#' Estimate a parsimonious leading-indicator model 
#' and evaluate short-horizon performance.
#' 
# ---- 6.1 Baseline: England ----
eng <- tsgc::england[, 1:2]

# Quick plot
mod2 <- tsgc::SSModelLeadingIndicator(eng, n.lag = 5)
p <- plot(
  mod2, title = "Daily COVID cases and Hospitalisations\n(England)",
  series.name.lead = "Cases", series.name.target = "Hospitalisations", take.log = TRUE
)
print(p)
# save_plot(p, "eng_hosp_lead_cases.png")

# Estimation window
est.start.eng <- as.Date("2021-04-30")
est.end.eng   <- as.Date("2021-07-24")
plt.len.eng   <- 14
n.lag         <- 4
n.forecasts   <- 7

# Define and estimate
out_eng <- tsgc::SSModelLeadingIndicator(
  Y = eng, n.lag = n.lag, q = NULL, LeadIndCol = 1, sea.period = 7,
  start.date = est.start.eng, end.date = est.end.eng
)
res_eng <- tsgc::estimate(out_eng)
summary(res_eng)

# Forecasts
p <- tsgc::plot_log_forecast(
  res_eng, Y = eng, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
  title = "Forecast of log growth rate of hospitalisation\n(England)"
); print(p)

p <- tsgc::plot_forecast(
  res_eng, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
  series.name = "hospitalisations", title = "Forecast of hospitalisation\n(England)"
); print(p)

p <- tsgc::plot_holdout(
  res_eng, Y = eng, n.ahead = n.forecasts,
  series.name = "hospitalisations", title = "Accuracy: Forecast of hospitalisation\n(England)"
); print(p)

if (SAVE_TABLES) {
  tsgc::write_results(
    res = res_eng, res.dir = tables_dir, n.ahead = n.forecasts,
    confidence.level = CONF_LEVEL, prefix = "england_hosp_lead_")
  message("Saved results for: england_hosp_lead")
}

#' ## Leading Indicator: Model with regressors
#' Augment with exogenous weather regressors for both lead and target, 
#' then forecast and evaluate.
#' 
## ---- 6.2 England With Regressors -xpred.  eval=TRUE ----
  xpred_lead <- xpred_targ <- england_weather_2021[, 1:4]
  mod_eng_x <- tsgc::SSModelLeadingIndicator(
    eng, n.lag = 4, xpred_lead = xpred_lead, xpred_targ = xpred_targ,
    start.date = est.start.eng, end.date = est.end.eng
  )
  res_eng_x <- tsgc::estimate(mod_eng_x)
  summary(res_eng_x)
  
  tsgc::supply_xpred.new(res_eng_x, england_weather_2021[, 1:4], idx = "lead")
  tsgc::supply_xpred.new(res_eng_x, england_weather_2021[, 1:4], idx = "targ")
  
  p <- tsgc::plot_log_forecast(
    res_eng_x, Y = eng, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
    title = "Forecast of log growth rate of hospitalisation\nwith regressors, England"
  ); print(p)
  
  p <- tsgc::plot_forecast(
    res_eng_x, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
    title = "Forecast of hospitalisation\nwith regressors, England", series.name = "hospitalisations"
  ); print(p)
  
  p <- tsgc::plot_holdout(
    res_eng_x, Y = eng, n.ahead = n.forecasts,
    title = "Accuracy: Forecast of hospitalisation\nwith regressors, England", series.name = "hospitalisations"
  ); print(p)
  
#' 
#' # Comparing Leading Indicator and Gompertz Models: UK & Italy
#' Compare UK forecasts from a UK-only Gompertz model 
#' vs a leading-indicator model with Italy as the lead across two windows. 
#' Both models are estimated on the same window with horizon 14 days 
#' and the uses same confidence level.
#' 
# ---- 7.1 UK vis-a-vis Italy ----
ukit <- tsgc::SSModelLeadingIndicator(tsgc::ukitaly, n.lag = 4)
p <- plot(
  ukit, title = "Daily COVID cases in UK and Italy",
  series.name.lead = "Italy", series.name.target = "UK", take.log = FALSE
)
print(p)
# save_plot(p, "ukit_cases_lead_cases.png")

# Case 1: First peak
n.forecasts <- 14; plt.length <- 30; CONF <- CONF_LEVEL
est.start <- as.Date("2020-02-25"); est.end <- as.Date("2020-04-01")
Yuk <- tsgc::ukitaly[, "UK"]

res_uk_gomp1 <- tsgc::estimate(
#tsgc::SSModelDynamicGompertz$new(Y = Yuk, q = q.default,
tsgc::SSModelDynamicGompertz(Y = Yuk, q = q.default,
                             start.date = est.start, 
                             end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_gomp1, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz)",
  plt.start.date = tail_date_minus(res_uk_gomp1$index, plt.length),
  series.name = "UK cases"
); print(p)

p <- tsgc::plot_holdout(
  res_uk_gomp1, Y = Yuk, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Accuracy: Forecast of daily COVID cases\nUK (Gompertz)", series.name = "UK cases"
); print(p)

# Leading indicator model
n.lag <- 14
res_uk_lead1 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(Y = tsgc::ukitaly, 
                                n.lag = n.lag, sea.period = 7,
                                start.date = est.start, 
                                end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_lead1, n.ahead = n.forecasts,
  title = "Leading indicator forecast\ndaily COVID cases in UK",
  plt.start.date = est.end - 30, series.name = "UK cases"
); print(p)

p <- tsgc::plot_holdout(
  res_uk_lead1, Y = tsgc::ukitaly, n.ahead = n.forecasts,
  title = "Accuracy: Leading indicator forecast\ndaily COVID cases in UK", series.name = "UK cases"
); print(p)

tsgc::plot_compare_forecast(list(res_uk_gomp1, res_uk_lead1), actual = tsgc::ukitaly[, "UK"])

## ---- Cross Validation ----
cv_models <- list()
# Model 1: Vanilla Gompertz
cv_models[["Vanilla_q"]]  <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, q = q.default, start.date = est.start, end.date = est.end
)
# Model 2: Vanilla Gompertz with AR(1)
cv_models[["Vanilla_ar1"]] <- tsgc::SSModelDynamicGompertz(
  Y = Yuk, start.date = est.start, end.date = est.end, ar1 = TRUE
)
# Model 3–6: Leading Indicator with different lags: 7, 10, 14, 18
for (i in c(7, 10, 14, 18)) {
  cv_models[[paste0("Lag", i)]] <- tsgc::SSModelLeadingIndicator(
    Y = tsgc::ukitaly, start.date = est.start, end.date = est.end, n.lag = i
  )
}
# Run cross-validation
tsgc::cross_val(
  Y = tsgc::ukitaly, model_list = cv_models, est.end.date = est.end,
  criterion="mape",
  n.estimate = 5, gap = 2
)

# Case 2: Extended window
est.start <- as.Date("2020-02-25"); est.end <- as.Date("2020-04-15")
res_uk_gomp2 <- tsgc::estimate(
#tsgc::SSModelDynamicGompertz$new(Y = Yuk, q = q.default,
tsgc::SSModelDynamicGompertz(Y = Yuk, q = q.default,
                             start.date = est.start, end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_gomp2, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz, extended)",
  plt.start.date = tail_date_minus(res_uk_gomp2$index, plt.length),
  series.name = "UK cases"
); print(p)

p <- tsgc::plot_holdout(
  res_uk_gomp2, Y = Yuk, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Accuracy: Forecast of daily COVID cases\nUK (Gompertz, extended)", series.name = "UK cases"
); print(p)

res_uk_lead2 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(Y = tsgc::ukitaly, n.lag = 14,
                                start.date = est.start, end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_lead2, n.ahead = n.forecasts,
  title = "Forecast of daily COVID cases\nUK (Leading indicator model, extended)",
  plt.start.date = est.end - plt.length, series.name = "UK cases"
); print(p)

p <- tsgc::plot_holdout(
  res_uk_lead2, Y = tsgc::ukitaly, n.ahead = n.forecasts,
  title = "Accuracy: Forecast of daily COVID cases\nUK (Leading indicator model, extended)", 
  series.name = "UK cases"
); print(p)

tsgc::plot_compare_forecast(list(res_uk_gomp2, res_uk_lead2), actual = tsgc::ukitaly[, "UK"])

#' 
#' # Extensions to Other Data Frequencies 
#' Demonstrate quarterly, monthly, and annual use-cases 
#' with appropriate seasonal settings and lead–lag structures.
#' 
#' ## Quarterly: Wii
#' Fit a quarterly Gompertz to Wii cumulative sales; forecast and evaluate.
#' 
## ---- 8.1 Quarterly: Wii ----
data(nintendo_sales, package = "tsgc")
wii <- nintendo_sales[, 1]

plot(wii)

# Note: The tsgc function requires the input cumulative series to be strictly 
# increasing in time. If the cumulative values exhibit plateaus—as in 
# the case of the Wii series—it is necessary to add small increments 
# to eliminate flat segments and allow model estimation. Note that a 
# strictly increasing cumulative series also implies that the underlying 
# (non-cumulative) series must be strictly positive.

# Model estimated for a strictly increasing segment
n.forecasts <- 4
est.start.q <- zoo::as.yearqtr("2006 Q4")
est.end.q   <- zoo::as.yearqtr("2010 Q3")

mod_wii <- tsgc::SSModelDynamicGompertz(  
  Y = wii, sea.period = 4, start.date = est.start.q, end.date = est.end.q
)
res_wii <- tsgc::estimate(mod_wii)
summary(res_wii)

# Cases with MA overlay
p <- plot(mod_wii, title = "Wii sales by quarter", 
          series.name = "Sales (Million)", MA_period = 4)
print(p)
# save_plot(p, "wii_sales_gomp_cases.png")

p <- tsgc::plot_log_forecast(res_wii, Y = wii, n.ahead = n.forecasts, title = "Log forecasts of Wii sales")
print(p)
# save_plot(p, "wii_sales_gomp_loggr_fcst.png")

p <- tsgc::plot_forecast(res_wii, n.ahead = n.forecasts, title = "Wii sales")
print(p)
# save_plot(p, "wii_sales_gomp_fcst.png")

p <- tsgc::plot_holdout(res_wii, Y = wii, n.ahead = n.forecasts, title = "Wii sales")
print(p)
# save_plot(p, "wii_sales_gomp_holdout.png")

# Explicit console+save example (reuse from Section 2 for consistency check)
p_holdout <- tsgc::plot_holdout(
  res_q, Y = cumulative_cases, n.ahead = n.forecasts, 
  confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases,\n14-days (Gauteng)", 
  series.name = "cases"
)
print(p_holdout)
# save_plot(p_holdout, "gauteng_cases_gomp_q005_holdout.png")

#' 
#' ## Leading indicator model with quarterly data: Switch to Wii
#' Model quarterly lead–lag between Switch and Wii; forecast and evaluate.
#' 
## ---- 8.2 Quarterly: Switch to Wii (Lead) ----
n.forecasts   <- 8
est.start.q2  <- zoo::as.yearqtr("2017 Q1")
est.end.q2    <- zoo::as.yearqtr("2019 Q4")
n.lag.q       <- zoo::as.yearqtr("2017 Q1") - zoo::as.yearqtr("2006 Q4")

y_q <- nintendo_sales[, c("wii", "switch_all")]
#mod_switch <- tsgc::SSModelLeadingIndicator$new(
mod_switch <- tsgc::SSModelLeadingIndicator(
  Y = y_q, sea.period = 4, n.lag = n.lag.q, start.date = est.start.q2, end.date = est.end.q2
)
res_switch <- tsgc::estimate(mod_switch)
summary(res_switch)

p <- tsgc::plot_log_forecast(res_switch, Y = y_q, n.ahead = n.forecasts, title = "Log forecasts of Switch sales")
print(p)
# save_plot(p, "switch_sales_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(res_switch, n.ahead = n.forecasts, title = "Switch sales", series.name = "sales")
print(p)
# save_plot(p, "switch_sales_lead_fcst.png")

p <- tsgc::plot_holdout(res_switch, Y = y_q, n.ahead = n.forecasts, title = "Switch sales", series.name = "sales")
print(p)
# save_plot(p, "switch_sales_lead_holdout.png")

#' 
#' ## Monthly: Plus500
#' Fit a monthly Gompertz to Plus500; forecast and evaluate.
#' 
## ---- 8.3 Monthly-etrading: Plus500 ----
data(etrading_apps, package = "tsgc")
Plus500 <- etrading_apps[, 1]

n.forecasts <- 4
est.start.m <- zoo::as.yearmon(2016)
est.end.m   <- zoo::as.yearmon(2021)

#mod_500 <- tsgc::SSModelDynamicGompertz$new(
mod_500 <- tsgc::SSModelDynamicGompertz(
  Y = Plus500, sea.period = 12, start.date = est.start.m, end.date = est.end.m
)
res_500 <- tsgc::estimate(mod_500)
summary(res_500)

p <- plot(mod_500, title = "Plus500 monthly downloads in France", series.name = "Monthly downloads", MA_period = 4)
print(p)
# save_plot(p, "plus500_downloads_gomp_cases.png")

p <- tsgc::plot_log_forecast(res_500, Y = Plus500, n.ahead = n.forecasts, title = "Log forecasts of Plus500 monthly downloads")
print(p)
# save_plot(p, "plus500_downloads_gomp_loggr_fcst.png")

p <- tsgc::plot_forecast(res_500, n.ahead = n.forecasts, title = "Plus500 monthly downloads")
print(p)
# save_plot(p, "plus500_downloads_gomp_fcst.png")

p <- tsgc::plot_holdout(res_500, Y = Plus500, n.ahead = n.forecasts, title = "Plus500 monthly downloads")
print(p)
# save_plot(p, "plus500_downloads_gomp_holdout.png")

#' 
#' ## Leading indicator model with momnthly data: DEGIRO to AvaTrade
#' Monthly leading-indicator model with DEGIRO as the lead for AvaTrade; 
#' forecast and evaluate.
#' 
## ---- 8.4 Monthly-leading: DEGIRO to AvaTrade (Lead) ----
n.forecasts <- 4
est.start.m2 <- zoo::as.yearmon(2017.5)
est.end.m2   <- zoo::as.yearmon(2021 + 1/12)
n.lag.m      <- zoo::as.yearmon(2017.5) - zoo::as.yearmon(2017)

y_m <- etrading_apps[, c("DEGIRO", "AvaTrade")]
#mod_500_lead <- tsgc::SSModelLeadingIndicator$new(
mod_500_lead <- tsgc::SSModelLeadingIndicator(
  Y = y_m, sea.period = 12, n.lag = n.lag.m, start.date = est.start.m2, end.date = est.end.m2
)
res_500_lead <- tsgc::estimate(mod_500_lead)
summary(res_500_lead)

p <- tsgc::plot_log_forecast(res_500_lead, Y = y_m, n.ahead = n.forecasts, title = "Log forecasts of AvaTrade monthly downloads")
print(p)
# save_plot(p, "avatrade_downloads_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(res_500_lead, n.ahead = n.forecasts, title = "AvaTrade monthly downloads", series.name = "downloads")
print(p)
# save_plot(p, "avatrade_downloads_lead_fcst.png")

p <- tsgc::plot_holdout(res_500_lead, Y = y_m, n.ahead = n.forecasts, title = "AvaTrade monthly downloads", series.name = "downloads")
print(p)
# save_plot(p, "avatrade_downloads_lead_holdout.png")

#' 
#' ## Annual: 3DS
#' Annual modelling with `sea.period = 0`; 
#' forecast and evaluate 3DS annual series.
#' 
## ---- 8.5 Annual: 3DS ----
n.forecasts <- 2
est.start.y <- zoo::as.yearmon(2011)
est.end.y   <- zoo::as.yearmon(2018)

# Convert quarterly to yearly (sample every 4th)
yearly_nintendo      <- nintendo_sales[4 * (1:19), c("wii", "3ds")]
threeds_xts          <- xts::xts(zoo::coredata(yearly_nintendo[, "3ds"]), order.by = zoo::yearmon(2005:2023))
yearly_nintendo_xts  <- xts::xts(zoo::coredata(yearly_nintendo), order.by = zoo::yearmon(2005:2023))

#mod_3ds <- tsgc::SSModelDynamicGompertz$new(
mod_3ds <- tsgc::SSModelDynamicGompertz(  
  Y = threeds_xts, sea.period = 0, start.date = est.start.y, end.date = est.end.y
)
res_3ds <- tsgc::estimate(mod_3ds)
summary(res_3ds)

p <- tsgc::plot_log_forecast(res_3ds, Y = threeds_xts, n.ahead = n.forecasts,
                             title = "Log forecasts of annual 3DS sales")
print(p)
# save_plot(p, "3ds_sales_gomp_loggr_fcst.png")

p <- tsgc::plot_forecast(res_3ds, n.ahead = n.forecasts, title = "Forecasts of annual 3DS sales")
print(p)
# save_plot(p, "3ds_sales_gomp_fcst.png")

p <- tsgc::plot_holdout(res_3ds, Y = threeds_xts, n.ahead = n.forecasts,
                        title = "Accuracy of annual 3DS sales forecasts")
print(p)
# save_plot(p, "3ds_sales_gomp_holdout.png")

#' 
#' ## Leading indicator model with annual data: Wii to 3DS (Lead)
#' Annual lead–lag model using Wii as lead for 3DS; forecast and evaluate.
#' 
## ---- 8.6 Annual-leading: Wii to 3DS (Lead) ----
n.lag.y <- zoo::as.yearmon(2011) - zoo::as.yearmon(2007)
mod_lead_y <- tsgc::SSModelLeadingIndicator(
  Y = yearly_nintendo_xts, sea.period = 0, n.lag = n.lag.y,
  start.date = est.start.y, end.date = est.end.y, LeadIndCol = 1
)
res_lead_y <- tsgc::estimate(mod_lead_y)
summary(res_lead_y)

p <- tsgc::plot_log_forecast(res_lead_y, Y = yearly_nintendo_xts, n.ahead = n.forecasts, title = "Log forecasts of 3DS annual sales")
print(p)
# save_plot(p, "3ds_sales_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(res_lead_y, n.ahead = n.forecasts, title = "Annual global 3DS sales", series.name = "sales (in Million)")
print(p)
# save_plot(p, "3ds_sales_lead_fcst.png")

p <- tsgc::plot_holdout(res_lead_y, Y = yearly_nintendo_xts, n.ahead = n.forecasts, title = "Annual global 3DS sales", series.name = "sales (in Million)")
print(p)
# save_plot(p, "3ds_sales_lead_holdout.png")

#' 
#' # Wrap-up
message("=== Run completed. Check 'results/Tables' and 'results/Images'. ===")
