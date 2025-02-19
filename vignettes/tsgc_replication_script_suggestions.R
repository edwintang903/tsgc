##############################################
# Replication Script: Time Series Growth Curves package `tsgc`
##############################################

# -----------------------------
# Library Management
# -----------------------------
library(tsgc)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggfortify)    
library(ggforce)
library(magrittr)
library(zoo)
library(latex2exp)
library(xts)
library(gridExtra)
library(here)
library(timetk)
library(ForecastComb) 
library(tidyr)

# -----------------------------
# Global Theme and Project Setup
# -----------------------------
# Set a consistent global theme.
theme_set(theme_economist_white(gray_bg = FALSE, base_size = 16))

# Define the base project path in one place. here() provides the root of your project.
base_path <- here() 

# -----------------------------
# Create Folder Structure with Error Handling
# -----------------------------
# Create a results folder with subfolders for Tables and Images.
results_dir <- file.path(base_path, 'results')
tables_dir  <- file.path(results_dir, 'Tables')
images_dir  <- file.path(results_dir, 'Images')

# Create the directories if they do not exist, with error handling.
for (dir_path in list(results_dir, tables_dir, images_dir)) {
  if (!dir.exists(dir_path)) {
    tryCatch(
      dir.create(dir_path, recursive = TRUE),
      error = function(e) {
        message("Error creating directory: ", dir_path, "\n", e)
      }
    )
  }
}

# -----------------------------
# 1. Gompertz Growth Curve Model
# 
# Parameter Definitions
# -----------------------------
# Gather all parameters in one centralised block for easy modification.
date.format      <- "%Y-%m-%d"
n.forecasts      <- 14
q                <- 0.005
confidence.level <- 0.68
plt.length       <- 30
estimation.date.start <- as.Date("2021-02-01")
estimation.date.end   <- as.Date("2021-04-19")

# -----------------------------
# Data Loading and Preparation: Gauteng Data
# -----------------------------
# Load Gauteng data (cumulative confirmed cases)
data(gauteng, package = "tsgc")
cumulative_cases <- gauteng[, 1]  

# -----------------------------
# Preliminary Examination of Data
# -----------------------------
# Calculate a centred 7-day moving average of daily differences.
ma.cent.new.cases <- zoo::rollmean(diff(cumulative_cases), 7, align = "center")
str(cumulative_cases)

# Identify the date with maximum new cases.
ma.cent.wave.3.idx.max <- tsgc::argmax(ma.cent.new.cases) %>% zoo::index()
ma.cent.wave.3.idx.max

# Prepare data for plotting by combining actual new cases and the moving average.
d <- cbind(diff(cumulative_cases), ma.cent.new.cases)
colnames(d) <- c('New Cases', 'Centered 7-day MA')
d.df <- data.frame(
  Date = index(d),
  New.Cases = coredata(d[, 1]),
  Centered.7.day.MA = coredata(d[, 2])
)

# -----------------------------
# Plotting: Daily Cases and Moving Average
# -----------------------------
# SUGGESTION: Consider creating a dedicated plotting function.
data_plot <- ggplot(data = d.df, aes(x = Date)) +
  geom_line(aes(y = New.Cases, color = "New Cases"), linewidth = 0.1) +
  geom_line(aes(y = Centered.7.day.MA, color = "Centered 7 day MA"), linewidth = 1) +
  scale_y_continuous(n.breaks = 10) +
  xlab("Day") +
  ylab("New cases") +
  scale_x_date(date_breaks = "60 days") +
  scale_color_manual(
    name = '',
    values = c('New Cases' = 'blue', 'Centered 7 day MA' = 'red')
  ) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.85),
    legend.title = element_text(size = 2),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold")
  )
data_plot

# -----------------------------
# Model Estimation Options for the Third Wave
# -----------------------------

idx.est <- (zoo::index(cumulative_cases) >= estimation.date.start) &
  (zoo::index(cumulative_cases) <= estimation.date.end)
y <- cumulative_cases[idx.est]

# -----------------------------
# Estimation: Diffuse Prior Model
# -----------------------------
# The signal-to-noise ratio was estimated as a free parameter in this step.
model_q <- SSModelDynamicGompertz$new(Y = y)
res_q <- estimate(model_q)
res_q


# -----------------------------
# Estimation: Fixed Signal-to-Noise Ratio Model
# -----------------------------
# Estimate model
model <- SSModelDynamicGompertz$new(Y = y, q = q)
res <- estimate(model)
res

# -----------------------------
# Forecasting: Log Growth Rate
# -----------------------------
tsgc::plot_log_forecast(
  res,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  plt.start.date = tail(res$index, 1) - plt.length,
  title = "Log Growth Rate Forecast"
)

# -----------------------------
# Forecasting: New Cases and Holdout Evaluation
# -----------------------------
tsgc::plot_new_cases(
  res,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  plt.start.date = tail(res$index, 1) - plt.length,
  series.name = "Cases"
)

tsgc::plot_holdout(
  res,
  Y = cumulative_cases,
  n.ahead = 14,
  confidence.level = confidence.level,
  date_format = date.format,
  series.name = "cases"
)

# Save the estimation results as CSV files.

tsgc::write_results(
  res = res,
  res.dir = results_dir,
  n.ahead = n.forecasts,
  confidence.level = confidence.level
)

# -----------------------------
# Reproduction Number Calculation
# -----------------------------
gen_int <- 4  # Generation interval in days

# Calculate reproduction number estimates and credible intervals.
r.t <- tail(exp(res$get_gy_ci() * gen_int), 7) %>% tk_tbl
r.t$name <- "Gauteng"
names(r.t) <- c("Date", "Rt", "lower", "upper", "name")
r.t

# Plot reproduction numbers.
res.rt <- ggplot(r.t, aes(x = Date)) +
  ylim(0, 1.4) +
  geom_line(aes(y = Rt, color = "Rt")) +
  geom_point(aes(y = Rt), color = "red", size = 3) +
  geom_segment(aes(xend = Date, yend = lower, y = Rt), color = "blue") +
  geom_segment(aes(xend = Date, yend = upper, y = Rt), color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "68%  Interval"), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "solid", linewidth = 1.5, color = "black") +
  scale_x_date(date_breaks = "1 day") +
  theme_light(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.2),
    legend.title = element_text(size = 2),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold")
  )
res.rt

# -----------------------------
# 2. Reinitialisation for a New Wave
#
# -----------------------------
# Update the estimation period.
estimation.date.end <- as.Date("2021-06-25")
idx.est <- (zoo::index(cumulative_cases) >= estimation.date.start) &
  (zoo::index(cumulative_cases) <= estimation.date.end)
y <- cumulative_cases[idx.est]

# Re-estimate the model over the new period.
model <- SSModelDynamicGompertz$new(Y = y, q = q)
res <- estimate(model)

# Trigger reinitialisation.
# Extract the smoothed slope and its standard deviation.
smoothed.slope.full <- xts::xts(res$output$alphahat[, "slope"], order.by = res$index)
smoothed.P.slope <- xts::xts(res$output$P[2, 2, -1], order.by = res$index)

# Calculate the smoothed slope and multiples of its standard error.
d2 <- cbind(smoothed.slope.full, sqrt(smoothed.P.slope),
            1.5 * sqrt(smoothed.P.slope),
            2 * sqrt(smoothed.P.slope))
d2.df <- data.frame(date = index(d2), coredata(d2))
colnames(d2.df) <- c("Date", "smoothed.slope", "sd.smoothed.slope",
                     "sd.smoothed.slope.1.5", "sd.smoothed.slope.2")

d2.df <- d2.df[d2.df$Date >= as.Date("2020-10-06"), ]

trigger.df <- d2.df %>%
  mutate(prev_smoothed.slope = lag(smoothed.slope)) %>%
  filter((smoothed.slope > sd.smoothed.slope.2 & prev_smoothed.slope < sd.smoothed.slope.2))
trigger.df

# Identify reinitialisation date based on a threshold condition.
reinit_zero.df <- d2.df %>%
  mutate(prev_smoothed.slope = lag(smoothed.slope)) %>%
  filter(Date < min(trigger.df$Date) & (smoothed.slope > 0 & prev_smoothed.slope < 0)) %>%
  arrange(desc(Date)) %>%
  slice(1)
reinit_zero.df

# Plot trigger and reinitialisation date.
ggplot(data = d2.df[d2.df$Date > '2021-02-11',], aes(x = Date)) +
  geom_line(aes(y = smoothed.slope, color = "smoothed.slope"), linewidth = 0.5) +
  geom_line(aes(y = sd.smoothed.slope, color = "sd.smoothed.slope"), linetype = "solid", linewidth = 0.25) +
  geom_line(aes(y = sd.smoothed.slope.1.5, color = "sd.smoothed.slope.1.5"), linetype = "solid", linewidth = 0.25) +
  geom_line(aes(y = sd.smoothed.slope.2, color = "sd.smoothed.slope.2"), linetype = "solid", linewidth = 0.5) +
  scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  geom_vline(data = trigger.df, aes(xintercept = Date), linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = reinit_zero.df, aes(xintercept = Date), linetype = "solid", size = 1, color = "black") +
  xlab("Day") +
  ylab("Slope") +
  scale_x_date(date_breaks = "10 days") +
  scale_color_manual(name = '',
                     values = c('smoothed.slope' = 'red',
                                'sd.smoothed.slope' = 'blue',
                                'sd.smoothed.slope.1.5' = 'green',
                                'sd.smoothed.slope.2' = 'black')) +
  theme_light(base_size = 11) +
  theme(
    legend.title = element_text(size = 2),
    legend.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(face = "bold")
  )

# Set the reinitialisation date.
reinit.dates <- "2021-04-21"

# Estimate the reinitialized model.
model <- SSModelDynamicGompertz$new(
  Y = y,
  q = q,
  reinit.date = as.Date(reinit.dates, format = date.format)
)
res.reinit <- estimate(model)

# -----------------------------
# Plotting: Forecasts after Reinitialisation
# -----------------------------
tsgc::plot_log_forecast(
  res.reinit,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  plt.start.date = tail(res.reinit$index, 1) - plt.length,
  title = "Forecast of ln(g_t) after reinitialisation."
)

tsgc::plot_new_cases(
  res.reinit,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  plt.start.date = tail(res.reinit$index, 1) - plt.length,
  series.name = "cases"
)

tsgc::plot_holdout(
  res.reinit,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  series.name = "cases"
)

# -----------------------------
# Holdout Evaluation Comparison
# -----------------------------
# Without reinitialisation.
tsgc::plot_holdout(
  res,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format
)

# -----------------------------
# 3. Leading Indicator Analysis: England Data
# -----------------------------

# Load England data and select the first two columns.
eng <- tsgc::england[, 1:2]

# Transform the data to calculate daily cases and log growth rates.
eng_full <- tsgc::add_daily_ldl(eng)
eng_daily <- eng_full[, 3:4]

# Plot daily new cases and admissions.
ggplot(log(eng_daily), aes(x = index(eng_daily))) +
  geom_line(aes(y = newCases, color = "Cases"), lwd = 0.85) +
  geom_line(aes(y = newAdmit, color = "Hospitalizations"), lwd = 0.85) +
  scale_color_manual(values = c("red", "blue"))+
  labs(
    title = "COVID Daily Cases and Hospitalizations in England",
    x = "Date",
    y = "log(Number)",
    color = "Legend"
  ) +
  theme_economist_white(gray_bg = FALSE, base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    text = element_text(size = rel(1.1)),
    axis.text = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1), margin = margin(r = 10)),
    axis.title.x = element_text(size = rel(1), margin = margin(t = 10)),
    plot.title = element_text(margin = margin(b = 5), face = "bold"),
    plot.caption = element_text(size = rel(1))
  ) +
  scale_x_date(labels = scales::date_format("%d %b %y"))


# Define estimation parameters for the leading indicator analysis.
estimation.date.start <- as.Date("2021-04-30")
estimation.date.end   <- as.Date("2021-07-24")
plt.length            <- 14  # Adjusted for this analysis
n.lag                 <- 4
n.forecasts           <- 7

# Select data for the estimation period.
idx.est <- (zoo::index(eng) >= estimation.date.start) &
  (zoo::index(eng) <= estimation.date.end)
y <- eng[idx.est, ]

# Recalculate daily data for the estimation sample.
eng_full <- tsgc::add_daily_ldl(y)
eng_daily <- eng_full[, 3:4]

# Plot daily new cases and admissions again.
ggplot(log(eng_daily), aes(x = index(eng_daily))) +
  geom_line(aes(y = newCases, color = "Cases"), lwd = 0.85) +
  geom_line(aes(y = newAdmit, color = "Hospitalizations"), lwd = 0.85) +
  scale_color_manual(values = c("red", "blue"))+
  labs(
    title = "COVID Daily Cases and Hospitalizations in England",
    x = "Date",
    y = "log(Number)",
    color = "Legend"
  ) +
  theme_economist_white(gray_bg = FALSE, base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    text = element_text(size = rel(1.1)),
    axis.text = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1), margin = margin(r = 10)),
    axis.title.x = element_text(size = rel(1), margin = margin(t = 10)),
    plot.title = element_text(margin = margin(b = 5), face = "bold"),
    plot.caption = element_text(size = rel(1))
  ) +
  scale_x_date(labels = scales::date_format("%d %b %y"))

# Define and estimate the leading indicator model.
out <- SSModelLeadingIndicator(Y = y, n.lag = n.lag, q = NA, LeadIndCol = 1, sea.period = 7)
res <- estimate(out)

# Plot forecasts.
plot_log_forecast(
  res,
  Y = eng,
  n.ahead = n.forecasts,
  plt.start.date = estimation.date.end - plt.length
)

plot_new_cases(
  res,
  n.ahead = n.forecasts,
  plt.start.date = estimation.date.end - plt.length,
  series.name = "hospitalizations"
)

plot_holdout(
  res,
  Y = eng,
  n.ahead = n.forecasts,
  series.name = "hospitalizations"
)

# Save results for the leading indicator analysis.
write_results(
  res = res,
  res.dir = results_dir,
  n.ahead = n.forecasts
)

# -----------------------------
# 4. Leading Indicator vs Gompertz Growth curves: UK and Italy Examples
# -----------------------------
# Example: UK-Italy Data analysis.
UKIT_full <- tsgc::add_daily_ldl(ukitaly)
ggplot(UKIT_full[, c(3, 4)], aes(x = index(UKIT_full))) +
  geom_line(aes(y = newCases, color = "Italy"), lwd = 0.85) +
  geom_line(aes(y = newAdmit, color = "UK"), lwd = 0.85) +
  scale_color_manual(values = c("red", "blue"))+
  labs(
    title = "COVID Daily Cases in UK and Italy",
    x = "Date",
    y = "Daily cases"
  ) +
  theme_economist_white(gray_bg = FALSE, base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    text = element_text(size = rel(1.1)),
    axis.text = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1), margin = margin(r = 10)),
    axis.title.x = element_text(size = rel(1), margin = margin(t = 10)),
    plot.title = element_text(margin = margin(b = 5), face = "bold"),
    plot.caption = element_text(size = rel(1))
  ) +
  scale_x_date(labels = scales::date_format("%d %b %y"))


# Case 1: First peak.
Y <- ukitaly[, "UK"]
estimation.date.start <- as.Date("2020-02-25")
n.forecasts <- 14
confidence.level <- 0.68
plt.length <- 30
estimation.date.end <- as.Date("2020-04-01")
idx.est <- (zoo::index(Y) >= estimation.date.start) &
  (zoo::index(Y) <= estimation.date.end)
y <- Y[idx.est]

model_q <- SSModelDynamicGompertz$new(Y = y, q = 0.005)
res <- estimate(model_q)

plot_new_cases(
  res,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  title = "UK predictions with vanilla growth model",
  date_format = date.format,
  plt.start.date = tail(res$index, 1) - plt.length,
  series.name = "UK cases"
)

plot_holdout(
  res,
  Y = Y[(tail(res$index, 1) + 0:n.forecasts)],
  title = "UK predictions with vanilla growth model",
  confidence.level = confidence.level,
  date_format = date.format,
  series.name = "UK cases"
)

# Compare to leading indicator model.
n.lag <- 14
n.forc <- 14
idx.est <- (zoo::index(ukitaly) >= estimation.date.start) &
  (zoo::index(ukitaly) <= estimation.date.end)
covid_xts <- ukitaly[idx.est]
out <- SSModelLeadingIndicator(Y = covid_xts, n.lag = n.lag, sea.period = 7)
res <- estimate(out)
plot_new_cases(
  res,
  n.ahead = n.forc,
  title = "UK predictions with leading indicator model",
  plt.start.date = estimation.date.end - 30,
  series.name = "UK cases"
)
plot_holdout(
  res,
  Y = ukitaly,
  title = "UK predictions with leading indicator model",
  n.ahead = n.forc,
  series.name = "UK cases"
)

# Case 2: Future peaks.
Y <- ukitaly[, "UK"]
estimation.date.start <- as.Date("2020-02-25")
n.forecasts <- 14
confidence.level <- 0.68
plt.length <- 30
estimation.date.end <- as.Date("2020-04-15")
idx.est <- (zoo::index(Y) >= estimation.date.start) &
  (zoo::index(Y) <= estimation.date.end)
y <- Y[idx.est]

model_q <- SSModelDynamicGompertz$new(Y = y, q = 0.005)
res <- estimate(model_q)

plot_new_cases(
  res,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  title = "UK predictions with vanilla growth model",
  date_format = date.format,
  plt.start.date = tail(res$index, 1) - plt.length,
  series.name = "UK cases"
)
plot_holdout(
  res,
  n.ahead = n.forecasts,
  Y = Y,
  title = "UK predictions with vanilla growth model",
  confidence.level = confidence.level,
  date_format = date.format,
  series.name = "UK cases"
)

# For leading indicator.
n.lag <- 14
n.forc <- 14
idx.est <- (zoo::index(ukitaly) >= estimation.date.start) &
  (zoo::index(ukitaly) <= estimation.date.end)
covid_xts <- ukitaly[idx.est]
out <- SSModelLeadingIndicator(Y = covid_xts, n.lag = n.lag)
res <- estimate(out)
plot_new_cases(
  res,
  n.ahead = n.forc,
  title = "UK predictions with leading indicator model",
  plt.start.date = estimation.date.end - plt.length,
  series.name = "UK cases"
)
plot_holdout(
  res,
  Y = ukitaly,
  title = "UK predictions with leading indicator model",
  n.ahead = n.forc,
  series.name = "UK cases"
)

# -----------------------------
# 5. Forecast Combination to Predict Hospitalisation from Multiple Lags of Cases
# -----------------------------
# Using the ForecastComb package, combine forecasts for hospitalisations.
Y <- england[, c("cum_cases", "cum_admissions")]
est.start.date <- as.Date("2020-09-01")
est.end.date <- as.Date("2020-10-30")
Y.reinit <- reinitialise_dataframe(Y, est.start.date)

#Split into 2 functions, rolling_forecasts and optimal_combine
comb_all<-combine_forecasts(
  Y.reinit,
  est.start.date,
  est.end.date,
  all_lags = c(2, 5, 7, 9),
  train_days = 20,
  test_days = 60,
  method = comb_BG,
  rolling = TRUE
)

est.start.date <- as.Date("2020-09-01")
est.end.date <- as.Date("2020-10-30")+80
idx.est <- (zoo::index(Y.reinit) >= est.start.date) &
  (zoo::index(Y.reinit) <= est.end.date)
y <- Y.reinit[idx.est]

future_preds<-matrix(nrow=14,ncol=4)
all_lags<-c(2, 5, 7, 9)
for (i in 1:4){
  j=all_lags[i]
  mod<-SSModelLeadingIndicator(Y=y, n.lag=j)
  resi<-estimate(mod)
  future_preds[,i]<-resi$predict_level(14)[,1]
}

### Are we writing the forecast combination results to file?
predict(comb_all,future_preds)