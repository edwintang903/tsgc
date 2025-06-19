##############################################
# Replication Script: Time Series Growth Curves package `tsgc`
##############################################

# -----------------------------
# Library Management
# -----------------------------
library(tsgc)
library(KFAS)
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
library(tidyr)
library(abind)

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
# Plotting: Daily Cases and Moving Average
# -----------------------------
# Get a glimpse of data by plotting its moving average series
mod1<-SSModelDynamicGompertz$new(Y=cumulative_cases)
plot(mod1, title="Gauteng daily cases", series.name="cases")

# -----------------------------
# Model Estimation Options for the Third Wave
# -----------------------------
y <- get_timeframe(cumulative_cases,estimation.date.start,estimation.date.end)

# -----------------------------
# Estimation: Diffuse Prior Model
# -----------------------------
# The signal-to-noise ratio was estimated as a free parameter in this step.
model <- SSModelDynamicGompertz$new(Y = y)
res <- estimate(model)
summary(res)

# -----------------------------
# Estimation: Diffuse Prior Model with AR(1) component
# -----------------------------
model_ar1 <- SSModelDynamicGompertz$new(Y = y, ar1=TRUE)
res_ar1 <- estimate(model_ar1)
summary(res_ar1)

# -----------------------------
# Estimation: Fixed Signal-to-Noise Ratio Model
# -----------------------------
# Estimate model
model_q <- SSModelDynamicGompertz$new(Y = y, q = q)
res_q <- estimate(model_q)
summary(res_q)

# -----------------------------
# Forecasting: Log Growth Rate
# -----------------------------
plot_log_forecast(res_q,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  plt.start.date = tail(res_q$index, 1) - plt.length,
  title = "Log Growth Rate Forecast of new cases (Gauteng)"
)

# -----------------------------
# Forecasting: New Cases and Holdout Evaluation
# -----------------------------
plot_new_cases(res_q,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  plt.start.date = tail(res_q$index, 1) - plt.length,
  title="14-day forecast for new cases (Gauteng)",
  series.name = "Cases"
)

plot_holdout(res_q,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  title="14-day forecast for new cases (Gauteng)",
  series.name = "cases"
)

# Save the estimation results as CSV files.

tsgc::write_results(
  res = res,
  res.dir = tables_dir,
  n.ahead = n.forecasts,
  confidence.level = confidence.level
)

# -----------------------------
# Estimation: Diffuse Prior Model with exogenous predictors
# -----------------------------
#Load Gauteng weather 
data(gauteng_weather_2021, package = "tsgc")
gauteng_weather<-gauteng_weather_2021[,c(1,3)]

# Set up model and estimate it
model_weather <- SSModelDynamicGompertz$new(Y = y, xpred=gauteng_weather)
res_weather <- estimate(model_weather)
summary(res_weather)

# Feed future weather data into the results object. Subsetting of gauteng_weather 
#is done inside the function.
supply_xpred.new(res_weather,gauteng_weather)

# Generate Forecasts
plot_log_forecast(res_weather,Y=cumulative_cases,n.ahead=n.forecasts,
                              plt.start.date = tail(res_weather$index, 1) - plt.length,
                              title = "Log Growth Rate Forecast")

plot_new_cases(res_weather,n.ahead=n.forecasts,
                           confidence.level = confidence.level,
                           date_format = date.format,
                           plt.start.date = tail(res_weather$index, 1) - plt.length,
                           title="14-day forecast for new cases Gauteng",
                           series.name = "Cases")

plot_holdout(res_weather,Y=cumulative_cases,n.ahead=n.forecasts,
                         confidence.level = confidence.level,
                         date_format = date.format,
                         title="14-day forecast for new cases Gauteng",
                         series.name = "cases")

plot_compare_forecast(list(res,res_q,res_ar1, res_weather), actual=cumulative_cases)

# -----------------------------
# Reproduction Number Calculation
# -----------------------------
gen_int <- 4  # Generation interval in days
ndays<-7 #Number of days to plot

# Calculate reproduction number estimates and credible intervals.
r.t <- estimate_r0(res_q, gen_int, ndays)
r.t

# Plot reproduction numbers.
estimate_r0(res_q, gen_int, ndays, show_plot = TRUE, 
            title="Gauteng Reproduction numbers")

# -----------------------------
# 2. Reinitialisation for a New Wave
#
# -----------------------------
# Update the estimation period.
estimation.date.end <- as.Date("2021-06-25")
y <- get_timeframe(cumulative_cases,estimation.date.start,estimation.date.end)

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
summary(res.reinit)

# # Estimate the reinitialized model with exogenous predictors.
# model.x <- SSModelDynamicGompertz$new(
#   Y = y,
#   xpred=gauteng_weather,
#   q = q,
#   reinit.date = as.Date(reinit.dates, format = date.format)
# )
# res.reinit.x <- estimate(model.x)
# summary(res.reinit.x)

# -----------------------------
# Plotting: Forecasts after Reinitialisation
# -----------------------------
plot_log_forecast(res.reinit,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  plt.start.date = tail(res.reinit$index, 1) - plt.length, 
  title = "Forecast of ln(g_t) after reinitialisation."
)

plot_new_cases(res.reinit,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  plt.start.date = tail(res.reinit$index, 1) - plt.length,
  title="With reinitialization",
  series.name = "cases"
)

plot_holdout(res.reinit,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  date_format = date.format,
  title="With reinitialization",
  series.name = "cases"
)

# -----------------------------
# Holdout Evaluation Comparison
# -----------------------------
# Without reinitialisation.
plot_holdout(res,
  Y = cumulative_cases,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  title="Without reinitialization",
  date_format = date.format
)

plot_compare_forecast(list(res,res.reinit), actual=cumulative_cases)

# -----------------------------
# 3. Leading Indicator Analysis: England Data
# -----------------------------
# Load England data and select the first two columns.
eng <- tsgc::england[, 1:2]

# Plot log daily new cases and admissions by calling the plot function.
mod2<-SSModelLeadingIndicator$new(eng, n.lag=5) #Choose any n.lag if only plotting is needed
plot(mod2,title="COVID Daily Cases and Hospitalizations in England",
          series.name.lead="Cases", 
          series.name.target="Hospitalizations",
          take.log=TRUE)

# Define estimation parameters for the leading indicator analysis.
estimation.date.start <- as.Date("2021-04-30")
estimation.date.end   <- as.Date("2021-07-24")
plt.length            <- 14  
n.lag                 <- 4
n.forecasts           <- 7

y <- get_timeframe(eng, estimation.date.start,estimation.date.end)

# Define the leading indicator model
out <- SSModelLeadingIndicator(Y = y, n.lag = n.lag, 
                               q = NULL, LeadIndCol = 1, sea.period = 7)

# Estimate the leading indicator model.
res <- estimate(out)
summary(res)

# Plot forecasts.
plot_log_forecast(res,
  Y = eng,
  n.ahead = n.forecasts,
  plt.start.date = estimation.date.end - plt.length,
  title="Forecasts of Log Growth rate of England hospitalizations"
)

plot_new_cases(res,
  n.ahead = n.forecasts,
  plt.start.date = estimation.date.end - plt.length,
  series.name = "hospitalizations",
  title="Forecasts of England hospitalizations"
)

plot_holdout(res,
  Y = eng,
  n.ahead = n.forecasts,
  series.name = "hospitalizations",
  title="Forecasts of England hospitalizations"
)

# Save results for the leading indicator analysis.
write_results(
  res = res,
  res.dir = results_dir,
  n.ahead = n.forecasts
)

# Cross-validation to identify the best n.lag
cross_val(y=eng[index(eng)>=estimation.date.start],
          est.end.date=estimation.date.end,n.ahead=7,all_lags=1:9,totaldays=3,
          vanilla=TRUE,freq=2,LeadIndCol=1, criterion="mape")

# -----------------------------
# Leading Indicator with exogenous predictors
# -----------------------------
xpred1<-xpred2<-england_weather_2021[,1:4]
mod<-SSModelLeadingIndicator$new(y, n.lag=4, xpred1=xpred1, xpred2=xpred2)
res_lead.x<-estimate(mod)
summary(res_lead.x)

supply_xpred.new(res_lead.x,england_weather_2021[,1:4],idx=1)
supply_xpred.new(res_lead.x,england_weather_2021[,1:4],idx=2)

# Plot forecasts
plot_new_cases(res_lead.x,
               n.ahead = n.forecasts,
               plt.start.date = estimation.date.end - plt.length,
               title="Forecasts of Log Growth rate of England hospitalizations")

plot_holdout(res_lead.x,
            Y = eng, n.ahead = n.forecasts,
            title="Forecasts of Log Growth rate of England hospitalizations"
)

plot_compare_forecast(res_lead.x)

# -----------------------------
# 4. Leading Indicator vs Gompertz Growth curves: UK and Italy Examples
# -----------------------------
# Example: UK-Italy Data analysis.
ukit<-SSModelLeadingIndicator$new(ukitaly, n.lag=4)
plot(ukit, title="COVID Daily Cases in UK and Italy",
          series.name.lead="Italy", 
          series.name.target="UK", take.log=FALSE)

# Case 1: First peak.
Y <- ukitaly[, "UK"]
estimation.date.start <- as.Date("2020-02-25")
n.forecasts <- 14
confidence.level <- 0.68
plt.length <- 30
estimation.date.end <- as.Date("2020-04-01")
y <- get_timeframe(Y, estimation.date.start, estimation.date.end)

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
covid_xts <- get_timeframe(ukitaly,estimation.date.start+1, estimation.date.end)
out <- SSModelLeadingIndicator(Y = covid_xts, n.lag = n.lag, sea.period = 7)
res_lead <- estimate(out)

plot_new_cases(res_lead,
  n.ahead = n.forc,
  title = "UK predictions with leading indicator model",
  plt.start.date = estimation.date.end - 30,
  series.name = "UK cases"
)

plot_holdout(res_lead,
  Y = ukitaly,
  title = "UK predictions with leading indicator model",
  n.ahead = n.forc,
  series.name = "UK cases"
)

plot_compare_forecast(list(res,res_lead), actual=ukitaly[,"UK"])

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
covid_xts <- get_timeframe(ukitaly,estimation.date.start,estimation.date.end)
out <- SSModelLeadingIndicator(Y = covid_xts, n.lag = n.lag)
res_lead <- estimate(out)
plot_new_cases(
  res_lead,
  n.ahead = n.forc,
  title = "UK predictions with leading indicator model",
  plt.start.date = estimation.date.end - plt.length,
  series.name = "UK cases"
)
plot_holdout(
  res_lead,
  Y = ukitaly,
  title = "UK predictions with leading indicator model",
  n.ahead = n.forc,
  series.name = "UK cases"
)

plot_compare_forecast(list(res,res_lead), actual=ukitaly[,"UK"])
