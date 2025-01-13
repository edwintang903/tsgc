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
#' @param res Results object estimated using the \code{estimate()} method.
#' @param n.ahead Number of forecasts (i.e. number of periods ahead to forecast
#' from end of estimation window). Default is 14.
#' @param confidence.level Width of prediction interval for \eqn{\ln g_t} to
#' use in forecasts of \eqn{y_t = \Delta Y_t}. Default is 0.68, which is
#' approximately one standard deviation for a Normal distribution.
#' @param date_format Date format. Default is \code{'\%Y-\%m-\%d'}.
#' @param title Title for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param plt.start.date First date of actual data (from estimation sample) to
#' plot on graph.\code{NULL} (i.e. plots all data in estimation window) by
#' default.
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
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#'
#' # Plot forecast of new cases 7 days ahead
#' plot_new_cases(res, n.ahead=7, confidence.level = 0.68, date_format = "%Y-%m-%d",
#' title = "Forecast new cases", plt.start.date = as.Date("2020-07-13"))
#'
#' @export
plot_new_cases <- function(object,...) {
  object$plot_new_cases(...)
}

#' @title Plots forecast and realised values of the log cumulative growth rate
#'
#' @description Plots actual and filtered values of the log cumulative growth
#' rate (\eqn{\ln(g_t)}) in the estimation sample and the forecast and realised
#' log cumulative growth rate out of the estimation sample.
#'
#' @param res Results object estimated using the \code{estimate()} method.
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
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' idx.eval <- (zoo::index(gauteng) >= as.Date("2020-07-20")) &
#'      zoo::index(gauteng) <= as.Date("2020-07-27")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#'
#' # Plot forecast and realised log growth rate of cumulative cases
#' plot_log_forecast(res, Y=gauteng, n.ahead = 7,
#'   title = "Forecast ln(g)", plt.start.date = as.Date("2020-07-13"))
#'
#' @export
plot_log_forecast <- function(object,...) {
  object$plot_log_forecast(...)
}

#' @title Plots the growth rates and slope of the log cumulative growth rate
#'
#' @description Plots the smoothed/filtered growth rate of the difference in
#' the cumulated variable (\eqn{g_y}), the smoothed/filtered growth rate of the
#' the cumulated variable (\eqn{g}), and the smoothed/filtered slope of
#' \eqn{\ln(g)}, \eqn{\gamma}.
#' Following Harvey and Kattuman (2021), we compute \eqn{g_{y,t}} as
#' \deqn{g_{y,t} = \exp(\delta_t) + \gamma_t.}
#'
#' @param res Results object estimated using the \code{estimate()} method.
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
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot filtered gy, g and gamma
#' plot_gy_components(res, plt.start.date = as.Date("2020-07-06"))
#'
#'@importFrom ggplot2 ggplot geom_line labs scale_x_date scale_y_continuous
#'@importFrom ggplot2 waiver
#'@importFrom ggplot2 theme margin scale_color_manual
#'@importFrom dplyr filter
#'@importFrom tidyr pivot_longer
#'@importFrom ggthemes theme_economist_white
#'@importFrom magrittr %>%
#'
#' @export
plot_gy_components <- function(object,...){
  object$plot_gy_components(...)
}

#' @title Plots the growth rates and slope of the log cumulative growth rate
#'
#' @description Plots the smoothed/filtered growth rate of the difference in the
#' cumulated variable (\eqn{g_y}) and the associated confidence intervals.
#'
#' @param res Results object estimated using the \code{estimate()} method.
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
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
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
plot_gy_ci <- function(object,...){
  object$plot_gy_ci(...)
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
#' @param res Results object estimated using the \code{estimate()} method.
#' @param Y Values of the cumulated variable, including the holdout sample.
#' sample (i.e. to which the forecasts should be compared to).
#' @param n.ahead The duration of the holdout sample.
#' @param confidence.level Width of prediction interval for \eqn{\ln(g_t)} to
#' use in forecasts of \eqn{y_t = \Delta Y_t}. Default is 0.68, which is
#' approximately one standard deviation for a Normal distribution.
#' @param series.name Name of the variable you are forecasting for the purposes
#' of a $y$-axis label. E.g. if \code{series.name = "Cases"} the \eqn{y}-axis
#' will show "New Cases".
#' @param date_format Date format, e.g. \code{'\%Y-\%m-\%d'}, which is the
#' default.
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
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' idx.eval <- (zoo::index(gauteng) >= as.Date("2020-07-20")) &
#'      zoo::index(gauteng) <= as.Date("2020-07-27")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#'
#' # Plot forecasts and outcomes over evaluation period
#' plot_holdout(object = res, Y = gauteng)
#'
#' @export
plot_holdout <- function(object,...) {
  object$plot_holdout(...)
}
