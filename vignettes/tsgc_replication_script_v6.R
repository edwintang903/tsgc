#############################################################
# Replication Script: Time Series Growth Curves `tsgc`
# Authors: Ashby, Harvey, Kattuman, Thamotheram, Tang 
# Last updated: <today>
#############################################################

## ================
## 0. Setup & Utils
## ================

# ---- Parameters & Toggles ----
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

# Estimation windows used in Section 1 by default
est.start.1 <- as.Date("2021-02-01")
est.end.1   <- as.Date("2021-04-19")

# ---- Quiet require + install ----
safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install it before running.", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- Load libraries  ----
libs <- c(
  "tsgc","KFAS","dplyr","ggplot2","ggthemes","ggfortify","ggforce",
  "magrittr","zoo","latex2exp","xts","gridExtra","here","timetk",
  "tidyr","abind","scales","grid"
)
invisible(lapply(libs, safe_library))

# ---- Global theme, seed, options ----
theme_set(ggthemes::theme_economist_white(gray_bg = FALSE, base_size = 16))
options(scipen = 7)  # fewer scientific notations

# ---- Project Paths  ----
setwd("C:/tsgc-incorporatedClasses")
# Ensure your working directory is the project root (where .Rproj or .here lives).
# base_path <- here::here()
base_path   <- getwd()
results_dir <- file.path(base_path, "results")
tables_dir  <- file.path(results_dir, "Tables")
images_dir  <- file.path(results_dir, "Images")

# ---- Directory creation ----
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- tryCatch({ dir.create(path, recursive = TRUE); TRUE },
                   error = function(e) { message("Failed to create: ", path, " :: ", e$message); FALSE })
    if (!ok) stop("Could not create required directory: ", path)
  }
}
invisible(lapply(list(results_dir, tables_dir, images_dir), ensure_dir))

# ---- Save/print helpers ----
safe_ggsave <- function(plot, filename, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI) {
  if (!SAVE_PLOTS) return(invisible(TRUE))
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved plot: ", normalizePath(filename, winslash = "/"))
}

# Always show on console; optionally save if fname supplied
save_plot <- function(p, fname = NULL) {
  print(p)
  if (!is.null(fname) && inherits(p, "ggplot")) {
    safe_ggsave(p, file.path(images_dir, fname))
  }
  invisible(p)
}

# Convenience for date windows used in plots
tail_date_minus <- function(index_vec, k) {
  if (length(index_vec) == 0) return(NA)
  tail(index_vec, 1) - k
}

## ============================
## 1. Vanilla Gompertz: Gauteng
## ============================

# Data
data(gauteng, package = "tsgc")
stopifnot(ncol(gauteng) >= 1)
cumulative_cases <- gauteng[, 1]

# Quick inspection (cases)
mod1 <- tsgc::SSModelDynamicGompertz$new(Y = cumulative_cases)
p_cases <- plot(mod1, title = "Gauteng daily cases", series.name = "cases")
save_plot(p_cases, "gauteng_cases_MA.png")

# ---- Estimation blocks ----
# 1a) Diffuse prior (no fixed q)
model_free <- tsgc::SSModelDynamicGompertz$new(
  Y = cumulative_cases,
  start.date = est.start.1,
  end.date   = est.end.1
)
res_free <- tsgc::estimate(model_free)
summary(res_free)

# 1b) Diffuse prior with AR(1)
model_ar1 <- tsgc::SSModelDynamicGompertz$new(
  Y = cumulative_cases, ar1 = TRUE,
  start.date = est.start.1,
  end.date   = est.end.1
)
res_ar1 <- tsgc::estimate(model_ar1)
summary(res_ar1)

# 1c) Fixed q
model_q <- tsgc::SSModelDynamicGompertz$new(
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1,
  end.date   = est.end.1
)
res_q <- tsgc::estimate(model_q)
summary(res_q)

# ---- Forecasts & evaluation (using res_q as reference) ----
n.forecasts <- n.forecasts.default
plt.length  <- plt.length.default

# Forecast of log growth rate 
p <- tsgc::plot_log_forecast(
  res_q, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_q$index, plt.length),
  title = "Forecast of log growth rate of cases (Gauteng)"
)
save_plot(p, "gauteng_cases_gomp_q005_loggr_fcst.png")

# New cases forecast
p <- tsgc::plot_forecast(
  res_q, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_q$index, plt.length),
  title = "Forecast of new cases\n14-days (Gauteng)", series.name = "Cases"
)
save_plot(p, "gauteng_cases_gomp_q005_fcst.png")

# Holdout evaluation
p <- tsgc::plot_holdout(
  res_q, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases\n14-days (Gauteng)",
  series.name = "cases"
)
save_plot(p, "gauteng_cases_gomp_q005_holdout.png")

if (SAVE_TABLES) {
  tsgc::write_results(res = res_free, res.dir = tables_dir, 
                      n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
                      prefix= "gauteng_gomp_q005_")
  message("Saved results for: gauteng_gomp_q005")
}

## ==================================================
## 1d. With Exogenous Predictors (Weather), Gauteng
## ==================================================
data(gauteng_weather_2021, package = "tsgc")
gauteng_weather <- gauteng_weather_2021[, c(1, 3)]

model_weather <- tsgc::SSModelDynamicGompertz$new(
  Y = cumulative_cases,
  xpred = gauteng_weather,
  start.date = est.start.1,
  end.date   = est.end.1
)
res_weather <- tsgc::estimate(model_weather)
summary(res_weather)

# Provide future xpred (has built-in checks)
tsgc::supply_xpred.new(res_weather, gauteng_weather)

# Forecasts (xpred)
p <- tsgc::plot_log_forecast(
  res_weather, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_weather$index, plt.length),
  title = "Forecast of log growth rate of cases\n(with regressors - weather)"
)
save_plot(p, "gauteng_cases_gomp_xpred_loggr_fcst.png")

p <- tsgc::plot_forecast(
  res_weather, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_weather$index, plt.length),
  title = "Forecast of new cases\nwith regressors (weather), Gauteng",
  series.name = "Cases"
)
save_plot(p, "gauteng_cases_gomp_xpred_fcst.png")

p <- tsgc::plot_holdout(
  res_weather, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases\nwith regressors", series.name = "cases"
)
save_plot(p, "gauteng_cases_gomp_xpred_holdout.png")

# Comparison (display only)
tsgc::plot_compare_forecast(list(res_free, res_q, res_ar1, res_weather), actual = cumulative_cases)

p <- tsgc::plot_compare_forecast(list(res_free, res_q, res_ar1, res_weather), actual = cumulative_cases)
print(p)

## ===========================
## Reproduction Number (R_t)
## ===========================
r.t <- tsgc::estimate_r0(res_q, gen_int, ndays)
if (SAVE_TABLES) {
  write.csv(r.t, row.names = FALSE, file = file.path(tables_dir, "gauteng_rt_gomp_q005.csv"))
  message("Saved gauteng_rt_gomp_q005.csv")
}

p <- tsgc::estimate_r0(res_q, gen_int, ndays, show_plot = TRUE, title = "Gauteng Reproduction numbers")
save_plot(p, "gauteng_rt_gomp_q005_plot.png")

## ==================================
## 2. Reinitialisation for a New Wave
## ==================================
est.end.2 <- as.Date("2021-06-25")

model_rei_base <- tsgc::SSModelDynamicGompertz$new(
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1, end.date = est.end.2
)
res_rei_base <- tsgc::estimate(model_rei_base)
summary(res_rei_base)

# KFS pieces via package accessors
kfs      <- tsgc::output(res_rei_base)
alphaHat <- tsgc::alphahat(kfs)   # T × m
Ptt_arr  <- tsgc::Ptt(kfs)        # m × m × T
idx      <- res_rei_base$index

# Smoothed slope (full length)
smthd.slpe.full <- xts::xts(alphaHat[, "slope"], order.by = idx)

# Variance of slope from Ptt (align lengths)
v_slope <- drop(Ptt_arr[2, 2, ])      # length T
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

# Plot trigger (show + save)
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
save_plot(p_trigger, "gauteng_cases_gomp_q005_reinit_trigger.png")

# Set reinit date (could also take from trigger.df/reinit_zero.df)
reinit.date <- as.Date("2021-04-21")

model_reinit <- tsgc::SSModelDynamicGompertz$new(
  Y = cumulative_cases, q = q.default,
  start.date = est.start.1, end.date = est.end.2,
  reinit.date = reinit.date
)
res_reinit <- tsgc::estimate(model_reinit)
summary(res_reinit)

# Forecasts after reinit
p <- tsgc::plot_log_forecast(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  plt.start.date = tail_date_minus(res_reinit$index, plt.length),
  title = "Forecast of log growth rate of cases\nafter reinitialisation"
)
save_plot(p, "gauteng_cases_gomp_q005_postreinit_loggr_fcst.png")

p <- tsgc::plot_forecast(
  res_reinit, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  plt.start.date = tail_date_minus(res_reinit$index, plt.length),
  title = "Forecast of new cases\nafter reinitialisation", series.name = "cases"
)
save_plot(p, "gauteng_cases_gomp_q005_postreinit_fcst.png")

p <- tsgc::plot_holdout(
  res_reinit, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL, 
  title = "Accuracy: Forecast of new cases\nwith reinitialization", series.name = "cases"
)
save_plot(p, "gauteng_cases_gomp_q005_postreinit_holdout.png")

# Compare holdout with/without reinitialization (baseline = res_rei_base)
tsgc::plot_holdout(
  res_rei_base, Y = cumulative_cases, n.ahead = n.forecasts,
  confidence.level = CONF_LEVEL, title = "Forecast without reinitialization (baseline window)"
)
tsgc::plot_compare_forecast(list(res_rei_base, res_reinit), actual = cumulative_cases)

## ==========================================
## 3. Leading Indicator: England (Daily)
## ==========================================
eng <- tsgc::england[, 1:2]

# Simple plot (cases)
mod2 <- tsgc::SSModelLeadingIndicator$new(eng, n.lag = 5)
p <- plot(
  mod2, title = "Daily COVID cases and Hospitalizations\n(England)",
  series.name.lead = "Cases", series.name.target = "Hospitalizations", take.log = TRUE
)
save_plot(p, "eng_hosp_lead_cases.png")

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
  title = "Forecast of log growth rate of hospitalization\n(England)"
)
save_plot(p, "eng_hosp_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(
  res_eng, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
  series.name = "hospitalizations", title = "Forecast of hospitalization\n(England)"
)
save_plot(p, "eng_hosp_lead_fcst.png")

p <- tsgc::plot_holdout(
  res_eng, Y = eng, n.ahead = n.forecasts,
  series.name = "hospitalizations", title = "Accuracy: Forecast of hospitalization\n(England)"
)
save_plot(p, "england_hosp_lead_holdout.png")

if (SAVE_TABLES) {
  tsgc::write_results(
    res = res_eng, res.dir = tables_dir, n.ahead = n.forecasts,
    confidence.level = CONF_LEVEL, prefix = "england_hosp_lead_")
  message("Saved results for: england_hosp_lead")
}

# ---- Leading indicator with regressors: xpred ----
if (exists("england_weather_2021")) {
  xpred1 <- xpred2 <- england_weather_2021[, 1:4]
  mod_eng_x <- tsgc::SSModelLeadingIndicator$new(
    eng, n.lag = 4, xpred1 = xpred1, xpred2 = xpred2,
    start.date = est.start.eng, end.date = est.end.eng
  )
  res_eng_x <- tsgc::estimate(mod_eng_x)
  summary(res_eng_x)
  
  tsgc::supply_xpred.new(res_eng_x, england_weather_2021[, 1:4], idx = 1)
  tsgc::supply_xpred.new(res_eng_x, england_weather_2021[, 1:4], idx = 2)
  
  p <- tsgc::plot_log_forecast(
    res_eng_x, Y = eng, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
    title = "Forecast of log growth rate of hospitalisation\nwith regressors, England"
  )
  save_plot(p, "england_hosp_lead_xpred_loggr_fcst.png")
  
  p <- tsgc::plot_forecast(
    res_eng_x, n.ahead = n.forecasts, plt.start.date = est.end.eng - plt.len.eng,
    title = "Forecast of hospitalisation\nwith regressors, England", series.name = "hospitalizations"
  )
  save_plot(p, "eng_hosp_lead_xpred_fcst.png")
  
  p <- tsgc::plot_holdout(
    res_eng_x, Y = eng, n.ahead = n.forecasts,
    title = "Accuracy: Forecast of hospitalization\nwith regressors, England", series.name = "hospitalizations"
  )
  save_plot(p, "england_hosp_lead_xpred_holdout.png")
} else {
  message("Skipping England xpred section: object 'england_weather_2021' not found.")
}

## =====================================================
## 4. Leading Indicator vs Gompertz: UK & Italy (Daily)
## =====================================================
ukit <- tsgc::SSModelLeadingIndicator$new(tsgc::ukitaly, n.lag = 4)
p <- plot(
  ukit, title = "Daily COVID cases in UK and Italy",
  series.name.lead = "Italy", series.name.target = "UK", take.log = FALSE
)
save_plot(p, "ukit_cases_lead_cases.png")

# Case 1: First peak
n.forecasts <- 14; plt.length <- 30; CONF <- 0.68
est.start <- as.Date("2020-02-25"); est.end <- as.Date("2020-04-01")
Yuk <- tsgc::ukitaly[, "UK"]

res_uk_gomp1 <- tsgc::estimate(
  tsgc::SSModelDynamicGompertz$new(Y = Yuk, q = q.default,
                                   start.date = est.start, end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_gomp1, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz)",
  plt.start.date = tail_date_minus(res_uk_gomp1$index, plt.length),
  series.name = "UK cases"
)
save_plot(p, "uk_cases_gomp_case1_fcst.png")

p <- tsgc::plot_holdout(
  res_uk_gomp1, Y = Yuk, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Accuracy: Forecast of daily COVID cases\nUK (Gompertz)", series.name = "UK cases"
)
save_plot(p, "uk_cases_gomp_case1_holdout.png")

# Leading indicator model
n.lag <- 14
res_uk_lead1 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(Y = tsgc::ukitaly, n.lag = n.lag, sea.period = 7,
                                start.date = est.start, end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_lead1, n.ahead = n.forecasts,
  title = "Leading indicator forecast\ndaily COVID cases in UK",
  plt.start.date = est.end - 30, series.name = "UK cases"
)
save_plot(p, "uk_cases_lead_case1_fcst.png")

p <- tsgc::plot_holdout(
  res_uk_lead1, Y = tsgc::ukitaly, n.ahead = n.forecasts,
  title = "Accuracy: Leading indicator forecast\ndaily COVID cases in UK", series.name = "UK cases"
)
save_plot(p, "uk_cases_lead_case1_holdout.png")

tsgc::plot_compare_forecast(list(res_uk_gomp1, res_uk_lead1), actual = tsgc::ukitaly[, "UK"])

# Case 2: Future peaks
est.start <- as.Date("2020-02-25"); est.end <- as.Date("2020-04-15")
res_uk_gomp2 <- tsgc::estimate(
  tsgc::SSModelDynamicGompertz$new(Y = Yuk, q = q.default,
                                   start.date = est.start, end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_gomp2, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Forecast of daily COVID cases\nUK (Gompertz, extended)",
  plt.start.date = tail_date_minus(res_uk_gomp2$index, plt.length),
  series.name = "UK cases"
)
save_plot(p, "uk_cases_gomp_case2_fcst.png")

p <- tsgc::plot_holdout(
  res_uk_gomp2, Y = Yuk, n.ahead = n.forecasts, confidence.level = CONF,
  title = "Accuracy: Forecast of daily COVID cases\nUK (Gompertz, extended)", series.name = "UK cases"
)
save_plot(p, "uk_cases_gomp_case2_holdout.png")

res_uk_lead2 <- tsgc::estimate(
  tsgc::SSModelLeadingIndicator(Y = tsgc::ukitaly, n.lag = 14,
                                start.date = est.start, end.date = est.end)
)

p <- tsgc::plot_forecast(
  res_uk_lead2, n.ahead = n.forecasts,
  title = "Forecast of daily COVID cases\nUK (Leading indicator model, extended)",
  plt.start.date = est.end - plt.length, series.name = "UK cases"
)
save_plot(p, "uk_cases_lead_case2_fcst.png")

p <- tsgc::plot_holdout(
  res_uk_lead2, Y = tsgc::ukitaly, n.ahead = n.forecasts,
  title = "Accuracy: Forecast of daily COVID cases\nUK (Leading indicator model, extended)", 
  series.name = "UK cases"
)
save_plot(p, "uk_cases_lead_case2_holdout.png")

tsgc::plot_compare_forecast(list(res_uk_gomp2, res_uk_lead2), actual = tsgc::ukitaly[, "UK"])

## ==================================
## 5. Other Data Resolutions
## ==================================

# -------- Quarterly: Nintendo (Wii) --------
data(nintendo_sales, package = "tsgc")
wii <- nintendo_sales[, 1]

n.forecasts <- 4
est.start.q <- zoo::as.yearqtr("2006 Q4")
est.end.q   <- zoo::as.yearqtr("2010 Q3")

mod_wii <- tsgc::SSModelDynamicGompertz$new(
  Y = wii, sea.period = 4, start.date = est.start.q, end.date = est.end.q
)
res_wii <- tsgc::estimate(mod_wii)
summary(res_wii)

# Cases with MA overlay
p <- plot(mod_wii, title = "Wii sales by quarter", series.name = "Sales (Million)", MA_period = 4)
save_plot(p, "wii_sales_gomp_cases.png")

p <- tsgc::plot_log_forecast(res_wii, Y = wii, n.ahead = n.forecasts, title = "Log forecasts of Wii sales")
save_plot(p, "wii_sales_gomp_loggr_fcst.png")

p <- tsgc::plot_forecast(res_wii, n.ahead = n.forecasts, title = "Wii sales")
save_plot(p, "wii_sales_gomp_fcst.png")

p <- tsgc::plot_holdout(res_wii, Y = wii, n.ahead = n.forecasts, title = "Wii sales")
save_plot(p, "wii_sales_gomp_holdout.png")

# Explicit console+save example
p_holdout <- tsgc::plot_holdout(
  res_q, Y = cumulative_cases, n.ahead = n.forecasts, confidence.level = CONF_LEVEL,
  title = "Accuracy: Forecast of new cases,\n14-days (Gauteng)", series.name = "cases"
)
save_plot(p_holdout, "gauteng_cases_gomp_q005_holdout.png")

# -------- Leading indicator: Switch vs Wii --------
n.forecasts   <- 8
est.start.q2  <- zoo::as.yearqtr("2017 Q1")
est.end.q2    <- zoo::as.yearqtr("2019 Q4")
n.lag.q       <- zoo::as.yearqtr("2017 Q1") - zoo::as.yearqtr("2006 Q4")

y_q <- nintendo_sales[, c("wii", "switch_all")]
mod_switch <- tsgc::SSModelLeadingIndicator$new(
  Y = y_q, sea.period = 4, n.lag = n.lag.q, start.date = est.start.q2, end.date = est.end.q2
)
res_switch <- tsgc::estimate(mod_switch)
summary(res_switch)

p <- tsgc::plot_log_forecast(res_switch, Y = y_q, n.ahead = n.forecasts, title = "Log forecasts of switch sales")
save_plot(p, "switch_sales_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(res_switch, n.ahead = n.forecasts, title = "Switch sales", series.name = "sales")
save_plot(p, "switch_sales_lead_fcst.png")

p <- tsgc::plot_holdout(res_switch, Y = y_q, n.ahead = n.forecasts, title = "Switch sales", series.name = "sales")
save_plot(p, "switch_sales_lead_holdout.png")

# -------- Monthly: eTrading Apps --------
data(etrading_apps, package = "tsgc")
Plus500 <- etrading_apps[, 1]

n.forecasts <- 4
est.start.m <- zoo::as.yearmon(2016)
est.end.m   <- zoo::as.yearmon(2021)

mod_500 <- tsgc::SSModelDynamicGompertz$new(
  Y = Plus500, sea.period = 12, start.date = est.start.m, end.date = est.end.m
)
res_500 <- tsgc::estimate(mod_500)
summary(res_500)

p <- plot(mod_500, title = "Plus500 monthly downloads in France", series.name = "Monthly downloads", MA_period = 4)
save_plot(p, "plus500_downloads_gomp_cases.png")

p <- tsgc::plot_log_forecast(res_500, Y = Plus500, n.ahead = n.forecasts, title = "Log forecasts of Plus500 monthly downloads")
save_plot(p, "plus500_downloads_gomp_loggr_fcst.png")

p <- tsgc::plot_forecast(res_500, n.ahead = n.forecasts, title = "Plus500 monthly downloads")
save_plot(p, "plus500_downloads_gomp_fcst.png")

p <- tsgc::plot_holdout(res_500, Y = Plus500, n.ahead = n.forecasts, title = "Plus500 monthly downloads")
save_plot(p, "plus500_downloads_gomp_holdout.png")

# -------- Leading indicator (monthly) --------
n.forecasts <- 4
est.start.m2 <- zoo::as.yearmon(2017.5)
est.end.m2   <- zoo::as.yearmon(2021 + 1/12)
n.lag.m      <- zoo::as.yearmon(2017.5) - zoo::as.yearmon(2017)

y_m <- etrading_apps[, c("DEGIRO", "AvaTrade")]
mod_500_lead <- tsgc::SSModelLeadingIndicator$new(
  Y = y_m, sea.period = 12, n.lag = n.lag.m, start.date = est.start.m2, end.date = est.end.m2
)
res_500_lead <- tsgc::estimate(mod_500_lead)
summary(res_500_lead)

p <- tsgc::plot_log_forecast(res_500_lead, Y = y_m, n.ahead = n.forecasts, title = "Log forecasts of AvaTrade monthly downloads")
save_plot(p, "avatrade_downloads_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(res_500_lead, n.ahead = n.forecasts, title = "AvaTrade monthly downloads", series.name = "downloads")
save_plot(p, "avatrade_downloads_lead_fcst.png")

p <- tsgc::plot_holdout(res_500_lead, Y = y_m, n.ahead = n.forecasts, title = "AvaTrade monthly downloads", series.name = "downloads")
save_plot(p, "avatrade_downloads_lead_holdout.png")

# -------- Yearly (via yearmon index): Nintendo 3DS --------
n.forecasts <- 2
est.start.y <- zoo::as.yearmon(2011)
est.end.y   <- zoo::as.yearmon(2018)

# Convert quarterly to yearly (sample every 4th)
yearly_nintendo      <- nintendo_sales[4 * (1:19), c("wii", "3ds")]
threeds_xts          <- xts::xts(zoo::coredata(yearly_nintendo[, "3ds"]), order.by = zoo::yearmon(2005:2023))
yearly_nintendo_xts  <- xts::xts(zoo::coredata(yearly_nintendo), order.by = zoo::yearmon(2005:2023))

mod_3ds <- tsgc::SSModelDynamicGompertz$new(
  Y = threeds_xts, sea.period = 0, start.date = est.start.y, end.date = est.end.y
)
res_3ds <- tsgc::estimate(mod_3ds)
summary(res_3ds)

p <- tsgc::plot_log_forecast(res_3ds, Y = threeds_xts, n.ahead = n.forecasts,
                             title = "Log Forecasts for upcoming annual EV sales in the US")
save_plot(p, "3ds_sales_gomp_loggr_fcst.png")

p <- tsgc::plot_forecast(res_3ds, n.ahead = n.forecasts, title = "Forecasts for upcoming annual EV sales in the US")
save_plot(p, "3ds_sales_gomp_fcst.png")

p <- tsgc::plot_holdout(res_3ds, Y = threeds_xts, n.ahead = n.forecasts,
                        title = "Accuracy of predictions for upcoming annual EV sales in the US")
save_plot(p, "3ds_sales_gomp_holdout.png")

# -------- Leading indicator (yearly) --------
n.lag.y <- zoo::as.yearmon(2011) - zoo::as.yearmon(2007)
mod_lead_y <- tsgc::SSModelLeadingIndicator$new(
  Y = yearly_nintendo_xts, sea.period = 0, n.lag = n.lag.y,
  start.date = est.start.y, end.date = est.end.y, LeadIndCol = 1
)
res_lead_y <- tsgc::estimate(mod_lead_y)
summary(res_lead_y)

p <- tsgc::plot_log_forecast(res_lead_y, Y = yearly_nintendo_xts, n.ahead = n.forecasts, title = "Log forecasts of 3ds sales")
save_plot(p, "3ds_sales_lead_loggr_fcst.png")

p <- tsgc::plot_forecast(res_lead_y, n.ahead = n.forecasts, title = "Annual global 3ds sales", series.name = "sales (in Million)")
save_plot(p, "3ds_sales_lead_fcst.png")

p <- tsgc::plot_holdout(res_lead_y, Y = yearly_nintendo_xts, n.ahead = n.forecasts, title = "Annual global 3ds sales", series.name = "sales (in Million)")
save_plot(p, "3ds_sales_lead_holdout.png")

message("=== Script completed. Check 'results/Tables' and 'results/Images'. ===")
