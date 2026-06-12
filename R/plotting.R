# Created by: Craig Thamotheram
# Created on: 11/02/2022

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


#' @title Plots the forecast of new cases (the difference of the cumulated
#' variable)
#'
#' @description Plots actual values of the difference in the cumulated variable,
#' the forecasts of the cumulated variable (both including and excluding the
#' seasonal component, where a seasonal is specified) and forecast intervals
#' around the forecasts. The forecast intervals are based on the prediction
#' intervals for \eqn{\ln(g_t)}.
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param n.ahead Number of forecasts (i.e. number of periods ahead to forecast
#' from end of estimation window). Default is 14.
#' @param confidence.level Width of prediction interval for \eqn{\ln g_t} to
#' use in forecasts of \eqn{y_t = \Delta Y_t}. Default is 0.68, which is
#' approximately one standard deviation for a Normal distribution.
#' @param title Title for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param plt.start.date First date of actual data (from estimation sample) to
#' plot on graph.\code{NULL} (i.e. plots all data in estimation window) by
#' default.
#' @param series.name The name of the series the growth rate is being computed
#' for. E.g. \code{'cases'}. Default is "target variable".
#'
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel scale_x_date
#' @importFrom ggplot2 geom_ribbon scale_size_manual margin
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom zoo coredata
#' @importFrom utils head tail
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#'
#' # Plot forecast of new cases 7 days ahead
#' plot_forecast(res, n.ahead=7, confidence.level = 0.68, 
#' title = "Forecast new cases", plt.start.date = as.Date("2020-07-13"),
#' series.name="cases")
#'
#' @export
plot_forecast <- function(res,n.ahead=7, confidence.level = 0.68, 
                           title=NULL, plt.start.date=NULL, 
                           series.name="target variable") {
  res$plot_forecast(n.ahead, confidence.level,
                        title, plt.start.date, 
                        series.name)
}

#' @title Plots forecast and realised values of the log cumulative growth rate
#'
#' @description Plots actual and filtered values of the log cumulative growth
#' rate (\eqn{\ln(g_t)}) in the estimation sample and the forecast and realised
#' log cumulative growth rate out of the estimation sample.
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param Y Cumulated dataset containing future values.
#' @param n.ahead The number of time periods ahead from the end of the sample
#' to be forecast. The default is 14.
#' @param plt.start.date Plot start date. Default is \code{NULL} which is the
#' start of the estimation sample.
#' @param title Plot title. Enter as text string.
#' @param caption Plot caption. Enter as text string.
#'
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel margin
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom utils tail
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' 
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#'
#' # Plot forecast and realised log growth rate of cumulative cases
#' plot_log_forecast(res, Y=gauteng, n.ahead = 7,
#'   title = "Forecast ln(g_t)", plt.start.date = as.Date("2020-07-13"))
#'
#' @export
plot_log_forecast <- function(res,Y, n.ahead = 14,
                              plt.start.date=NULL, title="", caption = "") {
  res$plot_log_forecast(Y, n.ahead,
                           plt.start.date, title, caption)
}

#' @title Plots the growth rates and slope of the log cumulative growth rate 
#' against the dates in estimation sample
#'
#' @description Plots the smoothed/filtered growth rate of the difference in
#' the cumulated variable (\eqn{g_y}), the smoothed/filtered growth rate of the
#' the cumulated variable (\eqn{g}), and the smoothed/filtered slope of
#' \eqn{\ln(g)}, \eqn{\gamma}.
#' Following Harvey and Kattuman (2021), we compute \eqn{g_{y,t}} as
#' \deqn{g_{y,t} = \exp(\delta_t) + \gamma_t.}
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param plt.start.date Plot start date. Default is \code{NULL} which is the
#' start of the estimation sample.
#' @param smoothed Logical value indicating whether to used the smoothed
#' estimates of \eqn{\delta} and \eqn{\gamma}. Default is \code{FALSE}, in
#' which case the filtered estimates are returned.
#' @param title Title for plot. Enter as text string. \code{NULL} (i.e. no
#' title) by default.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' 
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot filtered gy, g and gamma
#' plot_gy_components(res, plt.start.date = as.Date("2020-07-06"))
#'
#' @importFrom ggplot2 ggplot geom_line labs scale_x_date scale_y_continuous
#' @importFrom ggplot2 waiver
#' @importFrom ggplot2 theme margin scale_color_manual
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @importFrom ggthemes theme_economist_white
#' @importFrom magrittr %>%
#'
#' @export
plot_gy_components <- function(res,plt.start.date = NULL,
                               smoothed = FALSE, title = NULL){
  res$plot_gy_components(plt.start.date,smoothed, title)
}

#' @title Plots the growth rates and slope of the log cumulative growth rate
#'
#' @description Plots the smoothed/filtered growth rate of the difference in the
#' cumulated variable (\eqn{g_y}) and the associated confidence intervals.
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param plt.start.date Plot start date. Default is \code{NULL} which is the
#' start of the estimation sample.
#' @param smoothed Logical value indicating whether to used the smoothed
#' estimates of \eqn{\delta} and \eqn{\gamma}. Default is \code{FALSE}, in
#' which case the filtered estimates are returned.
#' @param title Title for plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param series.name The name of the series the growth rate is being computed
#' for. E.g. \code{'New cases'}.
#' @param pad.right Numerical value for the amount of time periods of blank
#' space you wish to leave on the right of the graph. Extends the horizontal
#' axis by the given number of time periods.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot filtered gy, g and gamma
#' plot_gy_ci(res, plt.start.date = as.Date("2020-07-13"))
#'
#' @importFrom ggplot2 ggplot geom_line geom_hline geom_ribbon labs
#' scale_color_manual scale_linetype_manual margin
#' @importFrom ggthemes theme_economist_white
#' @importFrom utils tail
#'
#' @export
plot_gy_ci <- function(res,plt.start.date = NULL, smoothed = FALSE,
                       title = NULL, series.name = NULL, pad.right = NULL){
  res$plot_gy_ci(plt.start.date, smoothed,
                    title, series.name, pad.right)
}

#' @title Plots the forecast of new cases (the difference of the cumulated
#' variable) over a holdout sample.
#'
#' @description Plots actual values of the difference in the cumulated variable,
#' the forecasts of the cumulated variable (both including and excluding the
#' seasonal component, where a seasonal is specified) and forecast intervals
#' around the forecasts, plus the actual outcomes from the holdout sample. The
#' forecast intervals are based on the prediction intervals for \eqn{\ln(g_t)}.
#' Also reports the mean absolute percentage prediction error over the holdout
#' sample.
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method. 
#' @param Y Values of the cumulated variable, including the holdout sample.
#' sample (i.e. to which the forecasts should be compared to).
#' @param n.ahead The duration of the holdout sample. Default is 14.
#' @param confidence.level Width of prediction interval for \eqn{\ln(g_t)} to
#' use in forecasts of \eqn{y_t = \Delta Y_t}. Default is 0.68, which is
#' approximately one standard deviation for a Normal distribution.
#' @param series.name Name of the variable you are forecasting for the purposes
#' of a $y$-axis label. E.g. if \code{series.name = "Cases"} the \eqn{y}-axis
#' will show "New Cases".
#' @param title Title for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param caption Caption for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no caption) by default.
#'
#' @importFrom xts as.xts
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel autoplot scale_x_date
#' @importFrom ggplot2 geom_ribbon scale_size_manual
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom zoo coredata index
#' @importFrom utils tail
#' @importFrom stats na.omit
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' n.ahead=7
#'
#' # Exapmle 1: Specify a Dynamic Gompertz model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot forecasts and outcomes over evaluation period
#' plot_holdout(res, Y = gauteng, n.ahead=n.ahead,series.name="cases")
#'
#' @export
plot_holdout <- function(res,Y, n.ahead=14,confidence.level = 0.68,
                         series.name = "target variable",
                         title= NULL, caption = NULL) {
  res$plot_holdout(Y, n.ahead, confidence.level, series.name, 
                      title, caption)
}

#' @title Forecast comparison plot
#'
#' @description Plots the forecasts of the cumulated variable (both including and excluding the
#' seasonal component, where a seasonal is specified) for selected models 
#' and optionally the actual outcomes from the holdout sample. The
#' forecast intervals are based on the prediction intervals for \eqn{\ln(g_t)}.
#'
#' @param results A list of `FilterResults` or `FilterResultsLI` objects, obtained from
#' \code{estimate()} method. 
#' @param actual Actual values of the cumulated variable, can be the raw dataset. 
#' Subsetting is done within the code.
#' @param sea.on Logical value indicating whether to plot the seasonality-adjusted
#' forecasts. Defaults to \code{TRUE}.
#' @param n.ahead The duration of the holdout sample. Default is 14.
#' @param title Title for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#'
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel autoplot scale_x_date
#' @importFrom ggplot2 geom_ribbon scale_size_manual
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom dplyr bind_rows
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' end.date<- as.Date("2020-07-20")
#' 
#' # Model 1: The signal-to-noise ratio is a free parameter
#' model <- SSModelDynamicGompertz$new(Y = gauteng, end.date=end.date)
#' res <- estimate(model)
#' 
#' # Model 2: The signal-to-noise ratio is fixed
#' model_q <- SSModelDynamicGompertz$new(Y = gauteng, end.date=end.date, q=0.005)
#' res_q <- estimate(model_q)
#' 
#' # Compare forecast
#' plot_compare_forecast(list(res,res_q), n.ahead=7)
#' 
#' @export
plot_compare_forecast <- function(results,  n.ahead = 14, sea.on = TRUE, actual = NULL,
                                  title = "Comparison of forecasts") {
  for (i in results){
   if (class(i)!="FilterResults" && class(i)!="FilterResultsLI"){
     stop("All elements in results list must be of the class FilterResults or FilterResultsLI.")
   }
  }
  
  # Automatically get object names as labels
  labels <- sapply(substitute(results)[-1], deparse)
  
  # Extract forecast data from each model
  prediction_list <- lapply(seq_along(results), function(i) {
    xts_pred <- results[[i]]$predict_level(n.ahead = n.ahead, sea.on = sea.on)[, 1]
    df_pred <- data.frame(
      date = as.Date(index(xts_pred)),
      forecast = as.numeric(xts_pred)
    )
    df_pred$model<- rep(labels[i], nrow(df_pred))
    return(df_pred)
  })
  
  # Combine all into one long-format data frame
  df_forecasts <- bind_rows(prediction_list)
  
  # Process actual values if provided
  if (!is.null(actual)) {
    actual<-na.omit(diff(actual))
    start.date<-tail(results[[1]]$index,1)+1
    end.date<-tail(results[[1]]$index,1)+n.ahead
    idx.est <- (zoo::index(actual) >= start.date) &
      (zoo::index(actual) <= end.date)
    actual<-actual[idx.est]
    
    if (inherits(actual, "xts")) {
      actual_df <- data.frame(
        date = as.Date(index(actual)),
        forecast = as.numeric(actual)
      )
    } else {
      actual_df <- data.frame(
        date = as.Date(names(actual)),
        forecast = as.numeric(actual)
      )
    }
    actual_df$model <- "Actual"
    df_forecasts <- dplyr::bind_rows(df_forecasts, actual_df)
  }
  df_forecasts$model<-unlist(df_forecasts$model)
  
  # Begin plot
  p <- ggplot(data = df_forecasts, aes(x = date)) +
    geom_line(aes(y = forecast, color = model),
              linewidth = 0.85) +
    labs(x = "Date", y = "Forecast", title = title) +
    theme_economist_white(gray_bg = FALSE, base_size = 14) +
    theme(legend.title = element_blank()) +
    theme(
      text = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      axis.title.y = element_text(size = rel(1), margin = margin(r = 10)),
      axis.title.x = element_text(size = rel(1), margin = margin(t = 10)),
      plot.title = element_text(margin = margin(b = 5)),
      plot.subtitle = element_text(size = rel(1), hjust = 0, margin = margin(t = 3))
    ) +
    scale_x_date(labels = scales::date_format("%d %b %y")) +
    scale_size_manual(values = c(1, 1.5, 1))
  
  return(p)
}
