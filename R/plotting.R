# Created by: Craig Thamotheram
# Created on: 11/02/2022
# Refactored: rewritten against idx_series/idx_calendar rather than
# xts/Date-indexed data. All plotting now lives in this single file as free
# functions (rather than scattered RefClass methods), operating on the
# idx_series/idx_calendar carried by FilterResults/FilterResultsLI/
# SSModelDynamicGompertz/SSModelLeadingIndicator objects.

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

## ----------------------------------------------------------------------
## Internal helpers: idx_series -> plottable data.frame
##
## All plot functions below go through these helpers so that there is a
## single place where the "position -> x-axis value" translation happens.
## If a calendar is available, the x-axis is real dates (via idx_to_date());
## otherwise it falls back to plain integer positions. Either way the
## resulting data.frame has an "x" column that ggplot's x aesthetic maps to,
## and an x-axis label/scale chosen to match.
## ----------------------------------------------------------------------

#' @keywords internal
#' @noRd
idx_series_df <- function(x, calendar = NULL) {
  stopifnot(is_idx_series(x))
  pos <- idx_positions(x)
  df <- as.data.frame(as.matrix(idx_values(x)))
  if (is_idx_calendar(calendar)) {
    df$x <- idx_to_date(calendar, pos)
  } else {
    df$x <- pos
  }
  df
}

#' @keywords internal
#' @noRd
idx_x_scale <- function(calendar = NULL) {
  if (is_idx_calendar(calendar) && isTRUE(calendar$posixct)) {
    ggplot2::scale_x_date(labels = scales::date_format("%d %b %y"))
  } else {
    ggplot2::scale_x_continuous()
  }
}

#' @keywords internal
#' @noRd
idx_x_lab <- function(calendar = NULL) {
  if (is_idx_calendar(calendar) && isTRUE(calendar$posixct)) "Date" else "Position"
}

## ----------------------------------------------------------------------
## SSModelDynamicGompertz / SSModelLeadingIndicator: plot(model)
## ----------------------------------------------------------------------

#' @title Plot the cumulated dataset underlying a fitted Gompertz model
#'
#' @description Plots the lagged differences of the cumulated dataset
#' \code{Y} stored in a \code{SSModelDynamicGompertz} object, optionally
#' overlaid with a centred moving average.
#'
#' @param x A \code{SSModelDynamicGompertz} object.
#' @param title Title for the plot. \code{NULL} (default) for no title.
#' @param series.name Name of the series being plotted, used in the y-axis
#' label. Default is \code{"target variable"}.
#' @param MA_period Number of positions in the centred moving average to
#' overlay. \code{0} or \code{1} disables the moving average. Default
#' \code{7}.
#' @param ... Unused.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs theme element_text
#' @importFrom ggplot2 scale_y_continuous
#'
#' @export
plot.SSModelDynamicGompertz <- function(x, title = NULL,
                                        series.name = "target variable",
                                        MA_period = 7, ...) {
  model <- x
  Y <- model$Y
  calendar <- model$calendar
  
  d <- idx_diff(Y, 1L)
  New.Cases <- NULL
  
  if (MA_period > 1) {
    vals <- idx_values(d)
    ma <- zoo::rollmean(vals, MA_period, align = "center")
    ma_pos <- idx_positions(d)[seq_along(ma) + (MA_period - 1) %/% 2]
    ma_series <- idx_series(ma, start = ma_pos[1])
    
    df <- idx_series_df(d, calendar)
    names(df)[1] <- "New.Cases"
    ma_df <- idx_series_df(ma_series, calendar)
    names(ma_df)[1] <- "Centered.MA"
    df$Centered.MA <- NA_real_
    df$Centered.MA[match(ma_df$x, df$x)] <- ma_df$Centered.MA
  } else {
    df <- idx_series_df(d, calendar)
    names(df)[1] <- "New.Cases"
  }
  
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = New.Cases, color = "New Cases"), linewidth = 0.1) +
    ggplot2::scale_y_continuous(n.breaks = 10) +
    ggplot2::labs(x = idx_x_lab(calendar), y = paste("New", series.name), title = title) +
    idx_x_scale(calendar) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 5),
      legend.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  if (MA_period > 1) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = Centered.MA, color = "Centered MA"), linewidth = 1) +
      ggplot2::scale_color_manual(name = '', values = c('Centered MA' = 'red'))
  }
  p
}

#' @title Plot the leading indicator and target variable underlying a fitted
#' Leading Indicator model
#'
#' @description Plots the daily incidence (new cases) of the leading
#' indicator and target variable stored in a \code{SSModelLeadingIndicator}
#' object, optionally on a log scale.
#'
#' @param x A \code{SSModelLeadingIndicator} object.
#' @param title Title for the plot. \code{NULL} (default) for no title.
#' @param series.name.lead Name of the leading indicator series, used in
#' the legend. Default \code{"Leading Indicator"}.
#' @param series.name.target Name of the target variable series, used in
#' the legend. Default \code{"Target Variable"}.
#' @param take.log Logical value indicating whether to plot the log of the
#' daily incidence. Default \code{TRUE}.
#' @param ... Unused.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs scale_color_manual theme
#' @importFrom ggplot2 element_text
#'
#' @export
plot.SSModelLeadingIndicator <- function(x, title = NULL,
                                         series.name.lead = "Leading Indicator",
                                         series.name.target = "Target Variable",
                                         take.log = TRUE, ...) {
  model <- x
  calendar <- model$calendar
  
  full <- add_daily_ldl(model$Y, LeadIndCol = model$LeadIndCol)
  newLead <- full$newLead
  newTarg <- full$newTarg
  
  common_pos <- intersect(idx_positions(newLead), idx_positions(newTarg))
  newLead <- newLead[common_pos]
  newTarg <- newTarg[common_pos]
  
  if (take.log) {
    newLead <- idx_series(log(idx_values(newLead)), start = newLead$start)
    newTarg <- idx_series(log(idx_values(newTarg)), start = newTarg$start)
    y.lab <- "log(Number)"
  } else {
    y.lab <- "Number"
  }
  
  df <- idx_series_df(newLead, calendar)
  names(df)[1] <- "newLead"
  df$newTarg <- idx_values(newTarg)
  
  ggplot2::ggplot(data = df, ggplot2::aes(x = x)) +
    ggplot2::labs(title = title, x = idx_x_lab(calendar), y = y.lab, color = "Legend") +
    ggplot2::geom_line(ggplot2::aes(y = newLead, color = series.name.lead), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = newTarg, color = series.name.target), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = c("red", "blue")) +
    idx_x_scale(calendar) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

## ----------------------------------------------------------------------
## plot_forecast
## ----------------------------------------------------------------------

#' @title Plots the forecast of new cases (the difference of the cumulated
#' variable)
#'
#' @description Plots actual values of the difference in the cumulated
#' variable, the forecast of the cumulated variable (including seasonal
#' components, where specified) and forecast intervals around the forecast.
#' The forecast intervals are based on the prediction intervals for
#' \eqn{\ln(g_t)}.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param n.ahead Number of forecasts (i.e. number of positions ahead to
#' forecast from the end of the estimation window). Default is 14.
#' @param confidence.level Width of the prediction interval for
#' \eqn{\ln g_t}, used in forecasts of \eqn{y_t = \Delta Y_t}. Default is
#' 0.68, approximately one standard deviation for a Normal distribution.
#' @param title Title for the plot. \code{NULL} (default) for no title.
#' @param plt.start First integer position of actual data (from the
#' estimation sample) to plot. \code{NULL} (default) plots all data in the
#' estimation window.
#' @param series.name Name of the series the growth rate is being computed
#' for, e.g. \code{'cases'}. Default \code{"target variable"}.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line geom_ribbon aes labs theme
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual
#' @importFrom ggplot2 element_blank element_text rel margin scale_size_manual
#'
#' @examples
#' library(tsgc)
#' set.seed(1)
#' Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1)
#' model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 100)
#' res <- estimate(model)
#' plot_forecast(res, n.ahead = 7, series.name = "cases")
#'
#' @export
plot_forecast <- function(res, n.ahead = 14, confidence.level = 0.68,
                          title = NULL, plt.start = NULL,
                          series.name = "target variable") {
  if (inherits(res, "FilterResultsLI")) {
    .plot_forecast_LI(res, n.ahead, confidence.level, title, plt.start, series.name)
  } else if (inherits(res, "FilterResults")) {
    .plot_forecast_gompertz(res, n.ahead, confidence.level, title, plt.start, series.name)
  } else {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
}

#' @keywords internal
#' @noRd
.plot_forecast_gompertz <- function(res, n.ahead, confidence.level, title,
                                    plt.start, series.name) {
  calendar <- res$calendar
  if (is.null(title)) title <- ""
  
  y.cum <- if (!is.null(res$reinit.idx)) {
    reinitialise_dataframe(res$data, res$reinit.idx)
  } else {
    res$data
  }
  if (is.null(plt.start)) plt.start <- idx_range(res$data)[1]
  
  y.hat.ci <- res$predict_level(n.ahead = n.ahead, confidence.level = confidence.level, sea.on = TRUE)
  
  d <- idx_diff(y.cum, 1L)
  keep <- idx_positions(d)[idx_positions(d) >= plt.start]
  d <- d[keep]
  
  df_plot <- idx_series_df(d, calendar)
  names(df_plot)[1] <- "Data"
  df_plot$Forecast <- NA_real_
  fc_df <- idx_series_df(idx_series(idx_values(y.hat.ci)[, 1], start = y.hat.ci$start), calendar)
  # Forecast positions may extend beyond the estimation data (they always
  # will, since this is a genuine out-of-sample forecast); assign forecast
  # values to matching existing rows, and append new rows for any forecast
  # positions not already present in df_plot.
  match_idx <- match(fc_df$x, df_plot$x)
  matched <- !is.na(match_idx)
  df_plot$Forecast[match_idx[matched]] <- fc_df[matched, 1]
  if (any(!matched)) {
    extra_rows <- data.frame(Data = NA_real_, x = fc_df$x[!matched], Forecast = fc_df[!matched, 1])
    df_plot <- rbind(df_plot, extra_rows[, names(df_plot)])
  }
  
  ci <- idx_series_df(idx_series(as.matrix(idx_values(y.hat.ci))[, 2:3, drop = FALSE], start = y.hat.ci$start), calendar)
  names(ci)[1:2] <- c("lower", "upper")
  
  Data <- Forecast <- lower <- upper <- NULL
  ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = Data, color = "Data"), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar), y = paste("New", series.name), title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar) +
    ggplot2::scale_size_manual(values = c(1, 1, 1))
}

#' @keywords internal
#' @noRd
.plot_forecast_LI <- function(res, n.ahead, confidence.level, title,
                              plt.start, series.name) {
  calendar <- res$calendar
  if (is.null(plt.start)) plt.start <- res$start
  
  sea <- res$predict_level(n.ahead = n.ahead, confidence.level = confidence.level, sea.on = TRUE)
  
  newTarg_full <- get_timeframe(res$data, res$start, res$end)
  newTarg_series <- idx_series(idx_values(newTarg_full)[, "newTarg"], start = newTarg_full$start)
  keep <- idx_positions(newTarg_series)[idx_positions(newTarg_series) >= plt.start]
  newTarg_series <- newTarg_series[keep]
  
  df_plot <- idx_series_df(newTarg_series, calendar)
  names(df_plot)[1] <- "newTarg"
  df_plot$Forecast <- NA_real_
  fc_df <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 1], start = sea$start), calendar)
  match_idx <- match(fc_df$x, df_plot$x)
  matched <- !is.na(match_idx)
  df_plot$Forecast[match_idx[matched]] <- fc_df[matched, 1]
  if (any(!matched)) {
    extra_rows <- data.frame(newTarg = NA_real_, x = fc_df$x[!matched], Forecast = fc_df[!matched, 1])
    df_plot <- rbind(df_plot, extra_rows[, names(df_plot)])
  }
  
  ci_plot <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 2:3, drop = FALSE], start = sea$start), calendar)
  names(ci_plot)[1:2] <- c("lwr", "upr")
  
  newTarg <- Forecast <- lwr <- upr <- NULL
  ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = newTarg, color = "Data"), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci_plot, ggplot2::aes(x = x, ymin = lwr, ymax = upr),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar), y = paste("New", series.name), title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar) +
    ggplot2::scale_size_manual(values = c(1, 1, 1))
}

#' @keywords internal
#' @noRd
idx_cbind_union <- function(...) {
  series_list <- list(...)
  all_pos <- sort(unique(unlist(lapply(series_list, idx_positions))))
  cols <- lapply(series_list, function(s) {
    vals <- idx_values(s)
    out <- rep(NA_real_, length(all_pos))
    out[match(idx_positions(s), all_pos)] <- vals
    out
  })
  mat <- do.call(cbind, cols)
  idx_series(mat, start = all_pos[1])
}

## ----------------------------------------------------------------------
## plot_log_forecast
## ----------------------------------------------------------------------

#' @title Plots forecast and realised values of the log cumulative growth
#' rate
#'
#' @description Plots actual and filtered values of the log cumulative
#' growth rate (\eqn{\ln(g_t)}) in the estimation sample and the forecast
#' and realised log cumulative growth rate out of the estimation sample.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param Y Cumulated dataset (\code{idx_series}) containing future values.
#' @param n.ahead Number of positions ahead from the end of the sample to
#' forecast. Default is 14.
#' @param plt.start Plot start position. \code{NULL} (default) is the
#' start of the estimation sample.
#' @param title Plot title.
#' @param caption Plot caption.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @export
plot_log_forecast <- function(res, Y, n.ahead = 14, plt.start = NULL,
                              title = "", caption = "") {
  if (inherits(res, "FilterResultsLI")) {
    .plot_log_forecast_LI(res, Y, n.ahead, plt.start, title, caption)
  } else if (inherits(res, "FilterResults")) {
    .plot_log_forecast_gompertz(res, Y, n.ahead, plt.start, title, caption)
  } else {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
}

#' @keywords internal
#' @noRd
.plot_log_forecast_gompertz <- function(res, Y, n.ahead, plt.start, title, caption) {
  calendar <- res$calendar
  model <- modelKFS(res$output)
  
  y.eval <- if (!is.null(res$reinit.idx)) {
    df2ldl(reinitialise_dataframe(Y, res$reinit.idx))
  } else {
    df2ldl(Y)
  }
  end.pos <- tail(res$index, 1)
  y.eval <- y.eval[idx_positions(y.eval)[idx_positions(y.eval) > end.pos]]
  
  y.series <- idx_series(as.numeric(gety(model)), start = res$index[1])
  p <- attr(model, "p")
  
  firstpred <- end.pos + 1L
  y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
  y.pred <- get_timeframe(y.hat.all$y.hat, firstpred, firstpred + n.ahead - 1L)
  filtered.level <- y.hat.all$level.t.t
  
  if (p != 1) stop("plot_log_forecast is only implemented for univariate (p=1) models.")
  
  actual_pos <- intersect(idx_positions(y.eval), idx_positions(y.pred))
  actual <- idx_series(idx_values(y.eval[actual_pos]), start = actual_pos[1])
  
  if (res$xpred_logical) {
    combined <- idx_cbind_union(y.series, idx_series(as.matrix(idx_values(y.pred))[, 1], start = y.pred$start), actual)
    colnames(combined$data) <- c("EstimationSample", "Forecast", "RealisedData")
    color_values <- c("Estimation\nSample" = 1, "Forecast" = 3, "Realised\nData" = "grey")
    linetype_values <- c("solid", "solid", "dashed")
  } else {
    combined <- idx_cbind_union(y.series, filtered.level,
                                idx_series(as.matrix(idx_values(y.pred))[, 1], start = y.pred$start), actual)
    colnames(combined$data) <- c("EstimationSample", "FilteredLevel", "Forecast", "RealisedData")
    color_values <- c("Estimation\nSample" = 1, "Filtered\nLevel" = 2, "Forecast" = 3, "Realised\nData" = "grey")
    linetype_values <- c("solid", "solid", "solid", "dashed")
  }
  
  if (!is.null(plt.start)) {
    keep <- idx_positions(combined)[idx_positions(combined) > plt.start]
    combined <- combined[keep]
  }
  
  df_plot <- idx_series_df(combined, calendar)
  
  EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
  p1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = EstimationSample, color = "Estimation\nSample"), linewidth = 0.85)
  
  if (!res$xpred_logical) {
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(y = FilteredLevel, color = "Filtered\nLevel"), linewidth = 0.85)
  }
  
  p1 +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = RealisedData, color = "Realised\nData"), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    idx_x_scale(calendar) +
    ggplot2::labs(x = idx_x_lab(calendar), y = "Log Growth Rate", caption = caption, title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    )
}

#' @keywords internal
#' @noRd
.plot_log_forecast_LI <- function(res, Y, n.ahead, plt.start, title, caption) {
  calendar <- res$calendar
  forcout_sea <- res$predict_all(n.ahead, sea.on = TRUE, return.all = FALSE)$y.hat
  old <- idx_series(idx_values(get_timeframe(res$data, res$start, res$end))[, "LDLtarg"],
                    start = res$start)
  
  full <- add_daily_ldl(Y, LeadIndCol = res$LeadIndCol)
  ldltarg_future <- full$LDLtarg
  ldltarg_future <- ldltarg_future[idx_positions(ldltarg_future)[idx_positions(ldltarg_future) > res$end]]
  n_take <- min(n.ahead, length(ldltarg_future))
  actual <- ldltarg_future[idx_positions(ldltarg_future)[seq_len(n_take)]]
  
  if (!any(res$xpred_logical)) {
    forcout <- res$predict_all(n.ahead, sea.on = FALSE, return.all = FALSE)$y.hat
    smldl_pred <- stats::predict(modelKFS(res$output), states = "trend")
    smldlh <- as.numeric(
      if (is.list(smldl_pred) && !is.data.frame(smldl_pred)) {
        smldl_pred[["LDLtarg"]]
      } else if (is.matrix(smldl_pred)) {
        smldl_pred[, "LDLtarg"]
      } else {
        smldl_pred
      }
    )
    filtered <- idx_series(smldlh, start = res$start)
    
    fc_series <- idx_series(as.matrix(idx_values(forcout))[, 1], start = forcout$start)
    forcout_full <- idx_rbind(filtered, fc_series[idx_positions(fc_series)[idx_positions(fc_series) > tail(idx_positions(filtered), 1)]])
    
    fc_sea_series <- idx_series(as.matrix(idx_values(forcout_sea))[, 1], start = forcout_sea$start)
    combined <- idx_cbind_union(old, forcout_full, fc_sea_series, actual)
    colnames(combined$data) <- c("EstimationSample", "FilteredLevel", "Forecast", "RealisedData")
    linetype_values <- c("solid", "solid", "solid", "dashed")
    color_values <- c(1, 2, 3, "grey")
  } else {
    fc_sea_series <- idx_series(as.matrix(idx_values(forcout_sea))[, 1], start = forcout_sea$start)
    combined <- idx_cbind_union(old, fc_sea_series, actual)
    colnames(combined$data) <- c("EstimationSample", "Forecast", "RealisedData")
    linetype_values <- c("solid", "solid", "dashed")
    color_values <- c(1, 3, "grey")
  }
  
  if (!is.null(plt.start)) {
    keep <- idx_positions(combined)[idx_positions(combined) >= plt.start]
    combined <- combined[keep]
  }
  df_plot <- idx_series_df(combined, calendar)
  
  EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
  p1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = EstimationSample, color = "Estimation\nSample"), linewidth = 0.85)
  
  if (!any(res$xpred_logical)) {
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(y = FilteredLevel, color = "Filtered\nLevel"), linewidth = 0.85)
  }
  
  p1 +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = RealisedData, color = "Realised\nData"), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    idx_x_scale(calendar) +
    ggplot2::labs(x = idx_x_lab(calendar), y = "Log Growth Rate", caption = caption, title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    )
}

## ----------------------------------------------------------------------
## plot_gy_components
## ----------------------------------------------------------------------

#' @title Plots the growth rate, level and slope of the log cumulative
#' growth rate
#'
#' @description Plots the smoothed/filtered growth rate of the difference
#' in the cumulated variable (\eqn{g_y}), the smoothed/filtered growth rate
#' of the cumulated variable (\eqn{g}), and the smoothed/filtered slope of
#' \eqn{\ln(g)}, \eqn{\gamma}. Following Harvey and Kattuman (2021), we
#' compute \eqn{g_{y,t}} as \deqn{g_{y,t} = \exp(\delta_t) + \gamma_t.}
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param plt.start Plot start position. \code{NULL} (default) is the
#' start of the estimation sample.
#' @param smoothed Logical value indicating whether to use the smoothed
#' estimates of \eqn{\delta} and \eqn{\gamma}. Default \code{FALSE}
#' (filtered estimates).
#' @param title Title for the plot. \code{NULL} (default) for no title.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line labs scale_x_date scale_y_continuous
#' @importFrom ggplot2 waiver theme margin scale_color_manual
#' @importFrom tidyr pivot_longer
#'
#' @export
plot_gy_components <- function(res, plt.start = NULL, smoothed = FALSE, title = NULL) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  calendar <- res$calendar
  if (is.null(plt.start)) plt.start <- res$index[1]
  
  gy.components <- res$get_growth_y(return.components = TRUE, smoothed = smoothed)
  gy.t <- gy.components[[1]]
  g.t <- gy.components[[2]]
  gamma.t <- gy.components[[3]]
  
  common_pos <- Reduce(intersect, lapply(list(gy.t, g.t, gamma.t), idx_positions))
  d <- idx_series(cbind(gy.t = idx_values(gy.t[common_pos]),
                        g.t = idx_values(g.t[common_pos]),
                        gamma.t = idx_values(gamma.t[common_pos])),
                  start = common_pos[1])
  
  keep <- idx_positions(d)[idx_positions(d) >= plt.start]
  d <- d[keep]
  
  df_plot <- idx_series_df(d, calendar)
  
  Variable <- Value <- NULL
  df_long <- tidyr::pivot_longer(df_plot, cols = c("gy.t", "g.t", "gamma.t"),
                                 names_to = "Variable", values_to = "Value")
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = Value, color = Variable)) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::facet_wrap(~ factor(Variable, c("gy.t", "g.t", "gamma.t")), ncol = 1, scales = "free_y") +
    ggplot2::labs(title = title, x = idx_x_lab(calendar), y = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#AA2045", "darkgrey", "black")) +
    idx_x_scale(calendar) +
    ggplot2::scale_y_continuous(breaks = ggplot2::waiver(), n.breaks = 4) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(b = 5)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      legend.position = "none"
    )
}

## ----------------------------------------------------------------------
## plot_gy_ci
## ----------------------------------------------------------------------

#' @title Plots the growth rate of the cumulative variable with confidence
#' intervals
#'
#' @description Plots the smoothed/filtered growth rate of the difference
#' in the cumulated variable (\eqn{g_y}) and the associated confidence
#' intervals.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param plt.start Plot start position. \code{NULL} (default) is the
#' start of the estimation sample.
#' @param smoothed Logical value indicating whether to use the smoothed
#' estimates. Default \code{FALSE} (filtered estimates).
#' @param title Title for the plot. \code{NULL} (default) for no title.
#' @param series.name Name of the series the growth rate is being computed
#' for, e.g. \code{'New cases'}.
#' @param pad.right Numerical value for the number of positions of blank
#' space to leave on the right of the graph.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line geom_hline geom_ribbon labs
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual margin
#'
#' @export
plot_gy_ci <- function(res, plt.start = NULL, smoothed = FALSE, title = NULL,
                       series.name = NULL, pad.right = NULL) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  calendar <- res$calendar
  if (is.null(plt.start)) plt.start <- res$index[1]
  
  gy.ci <- res$get_gy_ci(smoothed = smoothed)
  
  y.lab <- if (is.null(series.name)) "Growth rate" else paste("Growth rate of ", series.name, sep = "")
  
  df_plot <- idx_series_df(gy.ci, calendar)
  keep <- idx_positions(gy.ci) >= plt.start
  df_plot <- df_plot[keep, , drop = FALSE]
  
  fit <- upper <- lower <- NULL
  p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = fit), linewidth = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "green", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.3) +
    ggplot2::scale_color_manual(values = c("black")) +
    ggplot2::labs(title = title, x = idx_x_lab(calendar), y = y.lab) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    ) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "gray50", linewidth = 0.5)) +
    ggplot2::scale_linetype_manual(values = c("solid")) +
    idx_x_scale(calendar)
  
  if (!is.null(pad.right)) {
    if (is_idx_calendar(calendar) && isTRUE(calendar$posixct)) {
      end.date <- tail(df_plot$x, 1)
      p1 <- p1 + ggplot2::scale_x_date(limits = c(df_plot$x[1], end.date + pad.right))
    } else {
      end.pos <- tail(df_plot$x, 1)
      p1 <- p1 + ggplot2::scale_x_continuous(limits = c(df_plot$x[1], end.pos + pad.right))
    }
  }
  p1
}

## ----------------------------------------------------------------------
## plot_holdout
## ----------------------------------------------------------------------

#' @title Plots the forecast of new cases over a holdout sample
#'
#' @description Plots actual values of the difference in the cumulated
#' variable, the forecast of the cumulated variable (including seasonal
#' components, where specified) and forecast intervals around the
#' forecast, plus the actual outcomes from the holdout sample. Also reports
#' the mean absolute percentage error over the holdout sample.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param Y Values of the cumulated variable, including the holdout
#' sample.
#' @param n.ahead Duration of the holdout sample. Default is 14.
#' @param confidence.level Width of the prediction interval for
#' \eqn{\ln(g_t)}. Default 0.68.
#' @param series.name Name of the variable being forecast, for the y-axis
#' label.
#' @param title Title for the plot. \code{NULL} (default) for no title.
#' @param caption Caption for the plot. \code{NULL} (default) for no
#' caption.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line geom_ribbon labs theme
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual scale_size_manual
#'
#' @export
plot_holdout <- function(res, Y, n.ahead = 14, confidence.level = 0.68,
                         series.name = "target variable", title = NULL,
                         caption = NULL) {
  if (inherits(res, "FilterResultsLI")) {
    .plot_holdout_LI(res, Y, n.ahead, confidence.level, series.name, title, caption)
  } else if (inherits(res, "FilterResults")) {
    .plot_holdout_gompertz(res, Y, n.ahead, confidence.level, series.name, title, caption)
  } else {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
}

#' @keywords internal
#' @noRd
.plot_holdout_gompertz <- function(res, Y, n.ahead, confidence.level, series.name, title, caption) {
  calendar <- res$calendar
  estimation.end <- tail(res$index, 1)
  if (idx_range(Y)[2] <= estimation.end) {
    stop("Y must extend beyond the estimation end position to provide a holdout sample analysis.")
  }
  
  y.eval.diff <- idx_diff(Y, 1L)
  y.hat.ci <- res$predict_level(n.ahead = n.ahead, sea.on = TRUE, confidence.level = confidence.level)
  
  eval_pos <- idx_positions(y.eval.diff)
  keep_eval <- eval_pos[eval_pos > estimation.end & eval_pos < estimation.end + n.ahead + 1L]
  y.eval.diff <- y.eval.diff[keep_eval]
  
  common_pos <- intersect(idx_positions(y.eval.diff), idx_positions(y.hat.ci))
  d.eval <- data.frame(
    Actual = idx_values(y.eval.diff[common_pos]),
    Forecast = as.matrix(idx_values(y.hat.ci[common_pos]))[, 1]
  )
  d.eval <- stats::na.omit(d.eval)
  
  if (any(d.eval$Actual == 0)) {
    warning("Validation data contains zeros. MAPE is not a reliable measure.")
  }
  mape.sea <- signif(mean(100 * (abs(d.eval$Actual - d.eval$Forecast) / d.eval$Actual)), 4)
  smape <- round(mean(100 * (abs(d.eval$Actual - d.eval$Forecast) / (d.eval$Actual + d.eval$Forecast))), 2)
  mae <- signif(mean(abs(d.eval$Actual - d.eval$Forecast)), 4)
  rmse <- signif(sqrt(mean((d.eval$Actual - d.eval$Forecast)^2)), 4)
  
  df_plot <- idx_series_df(y.eval.diff, calendar)
  names(df_plot)[1] <- "Actual"
  fc_df <- idx_series_df(idx_series(as.matrix(idx_values(y.hat.ci))[, 1], start = y.hat.ci$start), calendar)
  df_plot$Forecast <- fc_df[match(df_plot$x, fc_df$x), 1]
  
  ci <- idx_series_df(idx_series(as.matrix(idx_values(y.hat.ci))[, 2:3, drop = FALSE], start = y.hat.ci$start), calendar)
  names(ci)[1:2] <- c("lower", "upper")
  
  Actual <- Forecast <- lower <- upper <- NULL
  ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = Actual, color = "Actual"), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar), y = paste("New", series.name), title = title,
                  subtitle = paste("MAPE: ", mape.sea, "%; SMAPE: ", smape, "%; MAE: ", mae, "; RMSE: ", rmse, ".", sep = "")) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(t = 3))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar) +
    ggplot2::scale_size_manual(values = c(1, 1.5, 1))
}

#' @keywords internal
#' @noRd
.plot_holdout_LI <- function(res, Y, n.ahead, confidence.level, series.name, title, caption) {
  calendar <- res$calendar
  if (idx_range(Y)[2] <= res$end) {
    stop("Y must extend beyond the estimation end position to provide a holdout sample analysis.")
  }
  
  sea <- res$predict_level(n.ahead = n.ahead, confidence.level = confidence.level, sea.on = TRUE)
  
  future_data <- add_daily_ldl(Y, LeadIndCol = res$LeadIndCol)
  newTarg_future <- future_data$newTarg
  newTarg_future <- newTarg_future[idx_positions(newTarg_future)[idx_positions(newTarg_future) > res$end]]
  if (n.ahead > length(newTarg_future)) {
    stop("The number of entries in the holdout sample is shorter than n.ahead. Please choose a smaller n.ahead, shorten the estimation period or provide more holdout data.")
  }
  actual <- newTarg_future[idx_positions(newTarg_future)[seq_len(n.ahead)]]
  
  common_pos <- intersect(idx_positions(actual), idx_positions(sea))
  compare <- data.frame(
    Actual = idx_values(actual[common_pos]),
    Forecast = as.matrix(idx_values(sea[common_pos]))[, 1]
  )
  
  if (any(compare$Actual == 0)) {
    warning("Validation data contains zeros. MAPE is not a reliable measure.")
  }
  mape.sea <- round(mean(100 * (abs(compare$Actual - compare$Forecast) / compare$Actual)), 2)
  smape <- round(mean(100 * (abs(compare$Actual - compare$Forecast) / (compare$Actual + compare$Forecast))), 2)
  mae <- signif(mean(abs(compare$Actual - compare$Forecast)), 3)
  rmse <- signif(sqrt(mean((compare$Actual - compare$Forecast)^2)), 3)
  
  df_plot <- idx_series_df(actual, calendar)
  names(df_plot)[1] <- "Actual"
  fc_df <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 1], start = sea$start), calendar)
  df_plot$Forecast <- fc_df[match(df_plot$x, fc_df$x), 1]
  
  ci_plot <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 2:3, drop = FALSE], start = sea$start), calendar)
  names(ci_plot)[1:2] <- c("lower", "upper")
  
  Actual <- Forecast <- lower <- upper <- NULL
  ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = Actual, color = "Actual"), linewidth = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci_plot, ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar), y = paste("New", series.name), title = title,
                  subtitle = paste("MAPE: ", mape.sea, "%; SMAPE: ", smape, "%; MAE: ", mae, "; RMSE: ", rmse, ".", sep = "")) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(t = 3))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar) +
    ggplot2::scale_size_manual(values = c(1, 1.5, 1))
}

## ----------------------------------------------------------------------
## plot_compare_forecast
## ----------------------------------------------------------------------

#' @title Forecast comparison plot
#'
#' @description Plots forecasts from a list of fitted models on the same
#' axes for visual comparison.
#'
#' @param results A list of \code{FilterResults} or \code{FilterResultsLI}
#' objects, obtained from the \code{estimate()} method.
#' @param n.ahead Duration of the holdout sample. Default is 14.
#' @param sea.on Logical value indicating whether to plot the
#' seasonality-adjusted forecasts. Defaults to \code{TRUE}.
#' @param actual An \code{idx_series} of actual cumulative values. The
#' function automatically extracts the observations corresponding to the
#' prediction period, so the complete series may be provided.
#' @param title Title for the plot. Defaults to \code{"Comparison of forecasts"}.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs theme element_blank
#' @importFrom ggplot2 element_text rel margin scale_size_manual
#'
#' @export
plot_compare_forecast <- function(results, n.ahead = 14, sea.on = TRUE,
                                  actual = NULL, title = "Comparison of forecasts") {
  for (r in results) {
    if (!inherits(r, "FilterResults") && !inherits(r, "FilterResultsLI")) {
      stop("All elements in results list must be of the class FilterResults or FilterResultsLI.")
    }
  }
  
  labels <- names(results)
  if (is.null(labels)) labels <- paste0("model", seq_along(results))
  
  calendar <- results[[1]]$calendar
  
  prediction_list <- lapply(seq_along(results), function(i) {
    pred <- results[[i]]$predict_level(n.ahead = n.ahead, sea.on = sea.on)
    fc <- idx_series(as.matrix(idx_values(pred))[, 1], start = pred$start)
    df <- idx_series_df(fc, calendar)
    names(df)[1] <- "forecast"
    df$model <- labels[i]
    df
  })
  df_forecasts <- do.call(rbind, prediction_list)
  
  if (!is.null(actual)) {
    actual.diff <- idx_diff(actual, 1L)
    end.pos <- tail(results[[1]]$index, 1)
    keep <- idx_positions(actual.diff)[idx_positions(actual.diff) > end.pos]
    keep <- keep[seq_len(min(n.ahead, length(keep)))]
    actual.diff <- actual.diff[keep]
    actual_df <- idx_series_df(actual.diff, calendar)
    names(actual_df)[1] <- "forecast"
    actual_df$model <- "Actual"
    df_forecasts <- rbind(df_forecasts, actual_df)
  }
  
  model <- forecast <- NULL
  ggplot2::ggplot(data = df_forecasts, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = forecast, color = model), linewidth = 0.85) +
    ggplot2::labs(x = idx_x_lab(calendar), y = "Forecast", title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(t = 3))
    ) +
    idx_x_scale(calendar) +
    ggplot2::scale_size_manual(values = c(1, 1.5, 1))
}