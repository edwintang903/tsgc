#' @title Options controlling how the x-axis of \code{idx_series} plots is
#' displayed
#'
#' @description Bundles the small set of cosmetic choices that decide how
#' plotting functions in this file turn integer \code{idx_series} positions
#' (and an optional \code{idx_calendar}) into an x-axis, and whether to
#' attach a small info box summarising the calendar. Pass an
#' \code{idx_axis_opts} object via the \code{axis} argument of the various
#' \code{plot.*}/\code{plot_*} functions; \code{NULL} (the default
#' everywhere) reproduces the previous behaviour exactly.
#'
#' @param mode One of:
#' \itemize{
#'   \item \code{"auto"} (default): \code{"date"} if a calendar with
#'   \code{posixct = TRUE} is supplied, else \code{"position"}.
#'   \item \code{"position"}: plain integer \code{idx_series} positions.
#'   \item \code{"steps"}: steps from the calendar's anchor position, i.e.
#'   \code{position - anchor_pos}. Labelled "Steps from \emph{anchor_name}" if
#'   the calendar has an \code{anchor_name}, else "Steps from anchor".
#'   \item \code{"time_since"}: the pattern-weighted calendar offset from
#'   the anchor (via \code{idx_calendar_offset()}), expressed in
#'   \code{calendar$unit}s. Labelled "Time since \emph{anchor_name} (unit)" or
#'   "Time since anchor (unit)". Also valid when \code{calendar$posixct}
#'   is \code{TRUE}: this simply plots the numeric offset instead of a real
#'   date, which can be more readable than dates for e.g. hourly/irregular
#'   series.
#'   \item \code{"date"}: real calendar dates via \code{idx_to_date()}
#'   (the original behaviour). Requires a calendar.
#' }
#' \code{"steps"}, \code{"time_since"} and \code{"date"} require a
#' calendar; they raise an error otherwise. \code{"position"} and
#' \code{"auto"} work with or without one.
#' @param info_box A single logical. If \code{TRUE}, adds a caption to the
#' plot summarising the calendar: anchor (and its name, if set), step size/
#' unit, and pattern. Defaults to \code{FALSE}. Ignored if no calendar is
#' supplied.
#' @param pattern_n \code{NULL} (default) to show the entire pattern in the
#' info box, or a single positive integer to truncate the displayed
#' pattern to its first \code{pattern_n} values (with a trailing "...").
#'
#' @returns An object of class \code{idx_axis_opts}.
#'
#' @examples
#' idx_axis_opts(mode = "steps", info_box = TRUE)
#' idx_axis_opts(mode = "time_since", info_box = TRUE, pattern_n = 5)
#'
#' @export
idx_axis_opts <- function(mode = c("auto", "position", "steps", "time_since", "date"),
                          info_box = FALSE, pattern_n = NULL) {
  mode <- match.arg(mode)
  if (length(info_box) != 1 || !is.logical(info_box) || is.na(info_box)) {
    stop("info_box must be a single logical (TRUE/FALSE).")
  }
  if (!is.null(pattern_n) &&
      (length(pattern_n) != 1 || !is.numeric(pattern_n) || pattern_n < 1)) {
    stop("pattern_n must be NULL or a single positive integer.")
  }
  structure(
    list(mode = mode, info_box = info_box,
         pattern_n = if (!is.null(pattern_n)) as.integer(pattern_n) else NULL),
    class = "idx_axis_opts"
  )
}

#' @keywords internal
#' @noRd
is_idx_axis_opts <- function(x) inherits(x, "idx_axis_opts")

#' @keywords internal
#' @noRd
idx_anchor_label <- function(calendar) {
  nm <- calendar$anchor_name
  if (!is.null(nm) && length(nm) == 1 && !is.na(nm) && nzchar(nm)) nm else "anchor"
}

#' @keywords internal
#' @noRd
idx_resolve_axis <- function(axis, calendar) {
  if (is.null(axis)) axis <- idx_axis_opts()
  stopifnot(is_idx_axis_opts(axis))
  mode <- axis$mode
  if (mode == "auto") {
    mode <- if (is_idx_calendar(calendar) && isTRUE(calendar$posixct)) "date" else "position"
  }
  if (mode %in% c("steps", "time_since", "date") && !is_idx_calendar(calendar)) {
    stop("axis mode '", mode, "' requires a calendar; none was supplied.")
  }
  if (mode == "time_since" && is_idx_calendar(calendar) && !is.null(calendar$multi_step)) {
    stop("axis mode 'time_since' is not supported for a calendar built via ",
         "idx_calendar_multi_step(): its steps cycle through heterogeneous ",
         "idx_step units that cannot be collapsed into a single numeric ",
         "offset. Use mode = 'steps' or 'date' instead.")
  }
  if (mode == "time_since" && is_idx_calendar(calendar) && is.na(calendar$amount)) {
    stop("axis mode 'time_since' is not supported for a calendar built via ",
         "idx_calendar_step(): it has no single 'amount' to express the ",
         "offset in. Use mode = 'steps' or 'date' instead.")
  }
  axis$mode <- mode
  axis
}

#' @keywords internal
#' @noRd
idx_series_df <- function(x, calendar = NULL, axis = NULL) {
  stopifnot(is_idx_series(x))
  pos <- idx_positions(x)
  df <- as.data.frame(as.matrix(idx_values(x)))
  axis <- idx_resolve_axis(axis, calendar)
  df$x <- switch(axis$mode,
                 date       = idx_to_date(calendar, pos),
                 steps      = pos - calendar$anchor_pos,
                 time_since = idx_calendar_offset(calendar, pos) * calendar$amount,
                 position   = pos
  )
  df
}

#' @keywords internal
#' @noRd
idx_x_scale <- function(calendar = NULL, axis = NULL) {
  axis <- idx_resolve_axis(axis, calendar)
  if (axis$mode == "date") {
    if (inherits(calendar$anchor, "POSIXct")) {
      ggplot2::scale_x_datetime(labels = scales::date_format("%d %b %y"),
                                breaks = scales::breaks_pretty(n = 8))
    } else if (inherits(calendar$anchor, "yearqtr")) {
      ggplot2::scale_x_continuous(
        breaks = scales::breaks_pretty(n = 8),
        labels = function(v) as.character(zoo::as.yearqtr(v))
      )
    } else if (inherits(calendar$anchor, "yearmon")) {
      ggplot2::scale_x_continuous(
        breaks = scales::breaks_pretty(n = 8),
        labels = function(v) as.character(zoo::as.yearmon(v))
      )
    } else {
      ggplot2::scale_x_date(labels = scales::date_format("%d %b %y"),
                            breaks = scales::breaks_pretty(n = 8))
    }
  } else {
    ggplot2::scale_x_continuous()
  }
}

#' @keywords internal
#' @noRd
idx_x_theme <- function(calendar = NULL, axis = NULL) {
  axis <- idx_resolve_axis(axis, calendar)
  if (axis$mode == "date") {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    ggplot2::theme()
  }
}

#' @keywords internal
#' @noRd
idx_forecast_x_theme <- function(calendar = NULL, axis = NULL) {
  axis <- idx_resolve_axis(axis, calendar)
  if (axis$mode == "date") {
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1,
        margin = ggplot2::margin(t = 6)
      ),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 18))
    )
  } else {
    ggplot2::theme()
  }
}

#' @keywords internal
#' @noRd
idx_x_lab <- function(calendar = NULL, axis = NULL) {
  axis <- idx_resolve_axis(axis, calendar)
  switch(axis$mode,
         date       = "Date",
         position   = "Position",
         steps      = paste0("Steps from ", idx_anchor_label(calendar)),
         time_since = paste0("Time since ", idx_anchor_label(calendar), " (", calendar$unit, ")")
  )
}

#' @keywords internal
#' @noRd
idx_info_box_caption <- function(calendar, axis) {
  if (is.null(axis) || !isTRUE(axis$info_box) || !is_idx_calendar(calendar)) return(NULL)
  
  anchor_line <- if (!is.null(calendar$anchor_name)) {
    paste0("Anchor: ", calendar$anchor_name, " (", format(calendar$anchor), ")")
  } else {
    paste0("Anchor: ", format(calendar$anchor))
  }
  
  if (!is.null(calendar$multi_step)) {
    step_lines <- vapply(seq_along(calendar$multi_step), function(i) {
      paste0("  [", i, "] ", idx_step_label(calendar$multi_step[[i]]))
    }, character(1))
    return(paste(
      anchor_line,
      paste0("Step: multi_step_pattern (", length(calendar$multi_step), " slots)"),
      paste(step_lines, collapse = "\n"),
      sep = "\n"
    ))
  }
  
  step_str <- if (is.na(calendar$amount)) {
    idx_step_label(calendar$step)
  } else {
    paste(calendar$amount, calendar$unit)
  }
  
  pat <- calendar$pattern
  pat_str <- if (!is.null(axis$pattern_n) && length(pat) > axis$pattern_n) {
    paste0(paste(pat[seq_len(axis$pattern_n)], collapse = ", "), ", ...")
  } else {
    paste(pat, collapse = ", ")
  }
  
  paste(
    anchor_line,
    paste0("Step: ", step_str),
    paste0("Pattern: ", pat_str),
    sep = "\n"
  )
}

#' @keywords internal
#' @noRd
idx_add_info_box <- function(p, calendar, axis) {
  cap <- idx_info_box_caption(calendar, axis)
  if (is.null(cap)) return(p)
  existing <- p$labels$caption
  full_cap <- if (!is.null(existing) && nzchar(existing)) paste(existing, cap, sep = "\n\n") else cap
  p + ggplot2::labs(caption = full_cap) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 7, family = "mono"))
}

## SSModelDynamicGompertz / SSModelLeadingIndicator: plot(model) ----

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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour (dates if a POSIXct calendar is present, else positions).
#' @param ... Unused.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs theme element_text
#' @importFrom ggplot2 scale_y_continuous scale_color_manual theme_bw element_blank
#'
#' @export
plot.SSModelDynamicGompertz <- function(x, title = NULL,
                                        series.name = "target variable",
                                        MA_period = 7, axis = NULL, ...) {
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
    
    df <- idx_series_df(d, calendar, axis)
    names(df)[1] <- "New.Cases"
    ma_df <- idx_series_df(ma_series, calendar, axis)
    names(ma_df)[1] <- "Centered.MA"
    df$Centered.MA <- NA_real_
    df$Centered.MA[match(ma_df$x, df$x)] <- ma_df$Centered.MA
  } else {
    df <- idx_series_df(d, calendar, axis)
    names(df)[1] <- "New.Cases"
  }
  
  p <- ggplot(data = df, aes(x = x)) +
    geom_line(aes(y = New.Cases), colour = "grey50", linewidth = 0.2, na.rm = TRUE) +
    scale_y_continuous(n.breaks = 10) +
    labs(x = idx_x_lab(calendar, axis), y = paste("New", series.name), title = title) +
    idx_x_scale(calendar, axis) +
    idx_x_theme(calendar, axis) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 5, l = 5)
    )
  if (MA_period > 1) {
    p <- p + geom_line(aes(y = Centered.MA, colour = "Centered MA"), linewidth = 1, na.rm = TRUE) + 
      scale_color_manual(values = c("Centered MA" = "red"))
  }
  idx_add_info_box(p, calendar, axis)
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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour.
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
                                         take.log = TRUE, axis = NULL, ...) {
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
  
  df <- idx_series_df(newLead, calendar, axis)
  names(df)[1] <- "newLead"
  df$newTarg <- idx_values(newTarg)
  
  p <- ggplot(data = df, aes(x = x)) +
    labs(title = title, x = idx_x_lab(calendar, axis), y = y.lab, color = "Legend") +
    geom_line(aes(y = newLead, color = series.name.lead), linewidth = 0.85, na.rm = TRUE) +
    geom_line(aes(y = newTarg, color = series.name.target), linewidth = 0.85, na.rm = TRUE) +
    scale_color_manual(values = c("red", "blue")) +
    idx_x_scale(calendar, axis) +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 5, r = 25, l = 5, b = 15)
    ) +
    idx_forecast_x_theme(calendar, axis)
  idx_add_info_box(p, calendar, axis)
}

## plot_forecast ----

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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour.
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
                          series.name = "target variable", axis = NULL) {
  if (inherits(res, "FilterResultsLI")) {
    .plot_forecast_LI(res, n.ahead, confidence.level, title, plt.start, series.name, axis)
  } else if (inherits(res, "FilterResults")) {
    .plot_forecast_gompertz(res, n.ahead, confidence.level, title, plt.start, series.name, axis)
  } else {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
}

#' @keywords internal
#' @noRd
.plot_forecast_gompertz <- function(res, n.ahead, confidence.level, title,
                                    plt.start, series.name, axis = NULL) {
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
  
  df_plot <- idx_series_df(d, calendar, axis)
  names(df_plot)[1] <- "Data"
  df_plot$Forecast <- NA_real_
  fc_df <- idx_series_df(idx_series(idx_values(y.hat.ci)[, 1], start = y.hat.ci$start), calendar, axis)
  match_idx <- match(fc_df$x, df_plot$x)
  matched <- !is.na(match_idx)
  df_plot$Forecast[match_idx[matched]] <- fc_df[matched, 1]
  if (any(!matched)) {
    extra_rows <- data.frame(Data = NA_real_, x = fc_df$x[!matched], Forecast = fc_df[!matched, 1])
    df_plot <- rbind(df_plot, extra_rows[, names(df_plot)])
  }
  
  ci <- idx_series_df(idx_series(as.matrix(idx_values(y.hat.ci))[, 2:3, drop = FALSE], start = y.hat.ci$start), calendar, axis)
  names(ci)[1:2] <- c("lower", "upper")
  
  Data <- Forecast <- lower <- upper <- NULL
  p <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = Data, color = "Data"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = paste("New", series.name), title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar, axis) +
    idx_forecast_x_theme(calendar, axis) +
    ggplot2::scale_size_manual(values = c(1, 1, 1))
  idx_add_info_box(p, calendar, axis)
}

#' @keywords internal
#' @noRd
.plot_forecast_LI <- function(res, n.ahead, confidence.level, title,
                              plt.start, series.name, axis = NULL) {
  calendar <- res$calendar
  if (is.null(plt.start)) plt.start <- res$start
  
  sea <- res$predict_level(n.ahead = n.ahead, confidence.level = confidence.level, sea.on = TRUE)
  
  newTarg_full <- get_timeframe(res$data, res$start, res$end)
  newTarg_series <- idx_series(idx_values(newTarg_full)[, "newTarg"], start = newTarg_full$start)
  keep <- idx_positions(newTarg_series)[idx_positions(newTarg_series) >= plt.start]
  newTarg_series <- newTarg_series[keep]
  
  df_plot <- idx_series_df(newTarg_series, calendar, axis)
  names(df_plot)[1] <- "newTarg"
  df_plot$Forecast <- NA_real_
  fc_df <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 1], start = sea$start), calendar, axis)
  match_idx <- match(fc_df$x, df_plot$x)
  matched <- !is.na(match_idx)
  df_plot$Forecast[match_idx[matched]] <- fc_df[matched, 1]
  if (any(!matched)) {
    extra_rows <- data.frame(newTarg = NA_real_, x = fc_df$x[!matched], Forecast = fc_df[!matched, 1])
    df_plot <- rbind(df_plot, extra_rows[, names(df_plot)])
  }
  
  ci_plot <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 2:3, drop = FALSE], start = sea$start), calendar, axis)
  names(ci_plot)[1:2] <- c("lwr", "upr")
  
  newTarg <- Forecast <- lwr <- upr <- NULL
  p <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = newTarg, color = "Data"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci_plot, ggplot2::aes(x = x, ymin = lwr, ymax = upr),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = paste("New", series.name), title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar, axis) +
    idx_forecast_x_theme(calendar, axis) +
    ggplot2::scale_size_manual(values = c(1, 1, 1))
  idx_add_info_box(p, calendar, axis)
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

## plot_log_forecast ----

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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour. If \code{info_box = TRUE}, the info box is appended below
#' any \code{caption} already supplied.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @export
plot_log_forecast <- function(res, Y, n.ahead = 14, plt.start = NULL,
                              title = "", caption = "", axis = NULL) {
  if (inherits(res, "FilterResultsLI")) {
    .plot_log_forecast_LI(res, Y, n.ahead, plt.start, title, caption, axis)
  } else if (inherits(res, "FilterResults")) {
    .plot_log_forecast_gompertz(res, Y, n.ahead, plt.start, title, caption, axis)
  } else {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
}

#' @keywords internal
#' @noRd
.plot_log_forecast_gompertz <- function(res, Y, n.ahead, plt.start, title, caption, axis = NULL) {
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
  
  df_plot <- idx_series_df(combined, calendar, axis)
  
  EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
  p1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = EstimationSample, color = "Estimation\nSample"), linewidth = 0.85, na.rm = TRUE)
  
  if (!res$xpred_logical) {
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(y = FilteredLevel, color = "Filtered\nLevel"), linewidth = 0.85, na.rm = TRUE)
  }
  
  p1 <- p1 +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = RealisedData, color = "Realised\nData"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    idx_x_scale(calendar, axis) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = "Log Growth Rate", caption = caption, title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    ) +
    idx_forecast_x_theme(calendar, axis)
  idx_add_info_box(p1, calendar, axis)
}

#' @keywords internal
#' @noRd
.plot_log_forecast_LI <- function(res, Y, n.ahead, plt.start, title, caption, axis = NULL) {
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
  df_plot <- idx_series_df(combined, calendar, axis)
  
  EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
  p1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = EstimationSample, color = "Estimation\nSample"), linewidth = 0.85, na.rm = TRUE)
  
  if (!any(res$xpred_logical)) {
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(y = FilteredLevel, color = "Filtered\nLevel"), linewidth = 0.85, na.rm = TRUE)
  }
  
  p1 <- p1 +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = RealisedData, color = "Realised\nData"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    idx_x_scale(calendar, axis) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = "Log Growth Rate", caption = caption, title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    ) +
    idx_forecast_x_theme(calendar, axis)
  idx_add_info_box(p1, calendar, axis)
}

## plot_gy_components ----

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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line labs scale_x_date scale_y_continuous
#' @importFrom ggplot2 waiver theme margin scale_color_manual
#' @importFrom tidyr pivot_longer
#'
#' @export
plot_gy_components <- function(res, plt.start = NULL, smoothed = FALSE, title = NULL, axis = NULL) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  calendar <- res$calendar
  if (is.null(plt.start)) plt.start <- idx_range(res$data)[1]
  
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
  
  df_plot <- idx_series_df(d, calendar, axis)
  
  Variable <- Value <- NULL
  df_long <- tidyr::pivot_longer(df_plot, cols = c("gy.t", "g.t", "gamma.t"),
                                 names_to = "Variable", values_to = "Value")
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = Value, color = Variable)) +
    ggplot2::geom_line(linewidth = 0.85, na.rm = TRUE) +
    ggplot2::facet_wrap(~ factor(Variable, c("gy.t", "g.t", "gamma.t")), ncol = 1, scales = "free_y") +
    ggplot2::labs(title = title, x = idx_x_lab(calendar, axis), y = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("#AA2045", "darkgrey", "black")) +
    idx_x_scale(calendar, axis) +
    idx_x_theme(calendar, axis) +
    ggplot2::scale_y_continuous(breaks = ggplot2::waiver(), n.breaks = 4) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(b = 5)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      legend.position = "none"
    )
  idx_add_info_box(p, calendar, axis)
}

## plot_gy_ci ----

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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line geom_hline geom_ribbon labs
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual margin
#'
#' @export
plot_gy_ci <- function(res, plt.start = NULL, smoothed = FALSE, title = NULL,
                       series.name = NULL, pad.right = NULL, axis = NULL) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  calendar <- res$calendar
  if (is.null(plt.start)) plt.start <- idx_range(res$data)[1]
  
  gy.ci <- res$get_gy_ci(smoothed = smoothed)
  
  y.lab <- if (is.null(series.name)) "Growth rate" else paste("Growth rate of ", series.name, sep = "")
  
  resolved_axis <- idx_resolve_axis(axis, calendar)
  df_plot <- idx_series_df(gy.ci, calendar, resolved_axis)
  keep <- idx_positions(gy.ci) >= plt.start
  df_plot <- df_plot[keep, , drop = FALSE]
  
  fit <- upper <- lower <- NULL
  p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = fit), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "green", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.3) +
    ggplot2::scale_color_manual(values = c("black")) +
    ggplot2::labs(title = title, x = idx_x_lab(calendar, resolved_axis), y = y.lab) +
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
    idx_x_scale(calendar, resolved_axis) +
    idx_x_theme(calendar, resolved_axis)
  
  if (!is.null(pad.right)) {
    if (resolved_axis$mode == "date") {
      end.date <- tail(df_plot$x, 1)
      if (inherits(calendar$anchor, "POSIXct")) {
        p1 <- p1 + ggplot2::scale_x_datetime(limits = c(df_plot$x[1], end.date + pad.right))
      } else {
        p1 <- p1 + ggplot2::scale_x_date(limits = c(df_plot$x[1], end.date + pad.right))
      }
    } else {
      end.pos <- tail(df_plot$x, 1)
      p1 <- p1 + ggplot2::scale_x_continuous(limits = c(df_plot$x[1], end.pos + pad.right))
    }
  }
  idx_add_info_box(p1, calendar, resolved_axis)
}

## plot_holdout ----

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
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line geom_ribbon labs theme
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual scale_size_manual
#'
#' @export
plot_holdout <- function(res, Y, n.ahead = 14, confidence.level = 0.68,
                         series.name = "target variable", title = NULL,
                         caption = NULL, axis = NULL) {
  if (inherits(res, "FilterResultsLI")) {
    .plot_holdout_LI(res, Y, n.ahead, confidence.level, series.name, title, caption, axis)
  } else if (inherits(res, "FilterResults")) {
    .plot_holdout_gompertz(res, Y, n.ahead, confidence.level, series.name, title, caption, axis)
  } else {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
}

#' @keywords internal
#' @noRd
.plot_holdout_gompertz <- function(res, Y, n.ahead, confidence.level, series.name, title, caption, axis = NULL) {
  calendar <- res$calendar
  if (!is_idx_series(Y)) {
    stop("Y must be an idx_series object.")
  }
  if (idx_ncol(Y) != 1) {
    stop("Y must be a single-column idx_series (subset a multi-column ",
         "series to the target column before passing it here).")
  }
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
    Actual = as.numeric(idx_values(y.eval.diff[common_pos])),
    Forecast = as.numeric(as.matrix(idx_values(y.hat.ci[common_pos]))[, 1])
  )
  d.eval <- stats::na.omit(d.eval)
  
  if (any(d.eval$Actual == 0)) {
    warning("Validation data contains zeros. MAPE is not a reliable measure.")
  }
  mape.sea <- signif(mean(100 * (abs(d.eval$Actual - d.eval$Forecast) / d.eval$Actual)), 4)
  smape <- round(mean(100 * (abs(d.eval$Actual - d.eval$Forecast) / (d.eval$Actual + d.eval$Forecast))), 2)
  mae <- signif(mean(abs(d.eval$Actual - d.eval$Forecast)), 4)
  rmse <- signif(sqrt(mean((d.eval$Actual - d.eval$Forecast)^2)), 4)
  
  df_plot <- idx_series_df(y.eval.diff, calendar, axis)
  names(df_plot)[1] <- "Actual"
  fc_df <- idx_series_df(idx_series(as.matrix(idx_values(y.hat.ci))[, 1], start = y.hat.ci$start), calendar, axis)
  df_plot$Forecast <- fc_df[match(df_plot$x, fc_df$x), 1]
  
  ci <- idx_series_df(idx_series(as.matrix(idx_values(y.hat.ci))[, 2:3, drop = FALSE], start = y.hat.ci$start), calendar, axis)
  names(ci)[1:2] <- c("lower", "upper")
  
  Actual <- Forecast <- lower <- upper <- NULL
  p <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = Actual, color = "Actual"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = paste("New", series.name), title = title,
                  caption = caption,
                  subtitle = paste("MAPE: ", mape.sea, "%; SMAPE: ", smape, "%; MAE: ", mae, "; RMSE: ", rmse, ".", sep = "")) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(t = 3)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar, axis) +
    idx_forecast_x_theme(calendar, axis) +
    ggplot2::scale_size_manual(values = c(1, 1.5, 1))
  idx_add_info_box(p, calendar, axis)
}

#' @keywords internal
#' @noRd
.plot_holdout_LI <- function(res, Y, n.ahead, confidence.level, series.name, title, caption, axis = NULL) {
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
    Actual = as.numeric(idx_values(actual[common_pos])),
    Forecast = as.numeric(as.matrix(idx_values(sea[common_pos]))[, 1])
  )
  
  if (any(compare$Actual == 0)) {
    warning("Validation data contains zeros. MAPE is not a reliable measure.")
  }
  mape.sea <- round(mean(100 * (abs(compare$Actual - compare$Forecast) / compare$Actual)), 2)
  smape <- round(mean(100 * (abs(compare$Actual - compare$Forecast) / (compare$Actual + compare$Forecast))), 2)
  mae <- signif(mean(abs(compare$Actual - compare$Forecast)), 3)
  rmse <- signif(sqrt(mean((compare$Actual - compare$Forecast)^2)), 3)
  
  df_plot <- idx_series_df(actual, calendar, axis)
  names(df_plot)[1] <- "Actual"
  fc_df <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 1], start = sea$start), calendar, axis)
  df_plot$Forecast <- fc_df[match(df_plot$x, fc_df$x), 1]
  
  ci_plot <- idx_series_df(idx_series(as.matrix(idx_values(sea))[, 2:3, drop = FALSE], start = sea$start), calendar, axis)
  names(ci_plot)[1:2] <- c("lower", "upper")
  
  Actual <- Forecast <- lower <- upper <- NULL
  p <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = Actual, color = "Actual"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Forecast, color = "Forecast"), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci_plot, ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = paste("New", series.name), title = title,
                  caption = caption,
                  subtitle = paste("MAPE: ", mape.sea, "%; SMAPE: ", smape, "%; MAE: ", mae, "; RMSE: ", rmse, ".", sep = "")) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(t = 3)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(1))
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
    idx_x_scale(calendar, axis) +
    idx_forecast_x_theme(calendar, axis) +
    ggplot2::scale_size_manual(values = c(1, 1.5, 1))
  idx_add_info_box(p, calendar, axis)
}

## plot_compare_forecast ----

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
#' @param labels Character vector of legend labels, one per element of
#' \code{results}, or \code{NULL} (default). If \code{NULL}, labels are
#' taken from \code{names(results)} if set (e.g.
#' \code{list(a = res1, b = res2)}), and otherwise default to the
#' deparsed expression used for each element in the \code{results} list
#' call (e.g. \code{plot_compare_forecast(list(res, res.reinit))} labels
#' the two lines \code{"res"} and \code{"res.reinit"}).
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour. Uses the calendar of \code{results[[1]]}.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs theme element_blank
#' @importFrom ggplot2 element_text rel margin scale_size_manual
#'
#' @export
plot_compare_forecast <- function(results, n.ahead = 14, sea.on = TRUE,
                                  actual = NULL, title = "Comparison of forecasts",
                                  labels = NULL, axis = NULL) {
  for (r in results) {
    if (!inherits(r, "FilterResults") && !inherits(r, "FilterResultsLI")) {
      stop("All elements in results list must be of the class FilterResults or FilterResultsLI.")
    }
  }
  
  if (!is.null(labels)) {
    if (length(labels) != length(results)) {
      stop("labels must have the same length as results.")
    }
  } else {
    labels <- names(results)
    if (is.null(labels)) {
      results_expr <- match.call()$results
      call_args <- if (is.call(results_expr)) as.list(results_expr)[-1] else NULL
      if (!is.null(call_args) && length(call_args) == length(results)) {
        arg_names <- names(call_args)
        labels <- vapply(seq_along(call_args), function(i) {
          nm <- arg_names[i]
          if (!is.null(nm) && nzchar(nm)) nm else deparse(call_args[[i]])
        }, character(1))
      }
    }
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- paste0("model", seq_along(results))
    }
  }
  
  calendar <- results[[1]]$calendar
  
  prediction_list <- lapply(seq_along(results), function(i) {
    pred <- results[[i]]$predict_level(n.ahead = n.ahead, sea.on = sea.on)
    fc <- idx_series(as.matrix(idx_values(pred))[, 1], start = pred$start)
    df <- idx_series_df(fc, calendar, axis)
    names(df)[1] <- "forecast"
    df$model <- labels[i]
    df
  })
  df_forecasts <- do.call(rbind, prediction_list)
  
  if (!is.null(actual)) {
    if (!is_idx_series(actual)) {
      stop("actual must be NULL or an idx_series object.")
    }
    if (idx_ncol(actual) != 1) {
      stop("actual must be a single-column idx_series (e.g. subset a ",
           "multi-column series to the target column before passing it ",
           "here); a multi-column actual produces a differently-shaped ",
           "data frame than the per-model forecasts and cannot be combined.")
    }
    actual.diff <- idx_diff(actual, 1L)
    end.pos <- tail(results[[1]]$index, 1)
    keep <- idx_positions(actual.diff)[idx_positions(actual.diff) > end.pos]
    keep <- keep[seq_len(min(n.ahead, length(keep)))]
    actual.diff <- actual.diff[keep]
    actual_df <- idx_series_df(actual.diff, calendar, axis)
    names(actual_df)[1] <- "forecast"
    actual_df$model <- "Actual"
    df_forecasts <- rbind(df_forecasts, actual_df)
  }
  
  model <- forecast <- NULL
  p <- ggplot2::ggplot(data = df_forecasts, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = forecast, color = model), linewidth = 0.85, na.rm = TRUE) +
    ggplot2::labs(x = idx_x_lab(calendar, axis), y = "Forecast", title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(1), margin = ggplot2::margin(t = 10)),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(t = 3))
    ) +
    idx_x_scale(calendar, axis) +
    idx_forecast_x_theme(calendar, axis) +
    ggplot2::scale_size_manual(values = c(1, 1.5, 1))
  idx_add_info_box(p, calendar, axis)
}

## plot_r0 ----

#' @title Plot the estimated reproduction number \eqn{R_t}
#'
#' @description Plots \eqn{R_t}, as computed by \code{\link{estimate_r0}},
#' together with its confidence interval, over the most recent
#' \code{n.ahead} integer positions of \code{res}. A horizontal reference
#' line is drawn at \eqn{R_t = 1}, the threshold separating a growing from
#' a shrinking epidemic/process.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param gen_int The mean generation interval, in the same integer-position
#' units as \code{res} (e.g. if positions are days, this is the mean
#' generation interval in days). Passed through to \code{estimate_r0()}.
#' @param n.ahead Number of most recent integer positions to plot, taken
#' from the end of \code{res}'s filtered/smoothed range. Default \code{7}.
#' Passed through to \code{estimate_r0()}.
#' @param smoothed Logical value indicating whether to use the smoothed
#' (\code{TRUE}) or filtered (\code{FALSE}, default) growth rate estimates
#' underlying \eqn{R_t}. Passed through to \code{estimate_r0()}.
#' @param confidence.level Confidence level for the confidence interval
#' around \eqn{R_t}. Default \code{0.68}. Passed through to
#' \code{estimate_r0()}.
#' @param title Title for the plot. \code{NULL} (default) for no title.
#' @param axis An \code{\link{idx_axis_opts}} object controlling x-axis
#' mode and info box, or \code{NULL} (default) for the previous default
#' behaviour (dates if a POSIXct calendar is present, else positions).
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' set.seed(1)
#' Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1)
#' model <- SSModelDynamicGompertz$new(Y = Y, q = NULL, end = 100)
#' res <- estimate(model)
#' plot_r0(res, gen_int = 5, n.ahead = 7)
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_segment geom_ribbon
#' @importFrom ggplot2 geom_hline labs theme scale_y_continuous
#'
#' @export
plot_r0 <- function(res, gen_int, n.ahead = 7, smoothed = FALSE,
                    confidence.level = 0.68, title = NULL, axis = NULL) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  
  calendar <- res$calendar
  r.t <- estimate_r0(res, gen_int = gen_int, n.ahead = n.ahead,
                     smoothed = smoothed, confidence.level = confidence.level)
  
  resolved_axis <- idx_resolve_axis(axis, calendar)
  r.t.series <- idx_series(
    as.matrix(r.t[, c("fit", "lower", "upper"), drop = FALSE]),
    start = r.t$Position[1]
  )
  df_plot <- idx_series_df(r.t.series, calendar, resolved_axis)
  
  fit <- lower <- upper <- NULL
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = fit, colour = "Rt"), na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = fit), colour = "red", size = 3) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = lower, y = fit), colour = "blue") +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = upper, y = fit), colour = "blue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper,
                                      fill = paste0(round(confidence.level * 100), "% Interval")),
                         alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "solid", linewidth = 1, colour = "black") +
    ggplot2::labs(title = title, x = idx_x_lab(calendar, resolved_axis),
                  y = expression(R[t])) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::theme_light(base_size = 12) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.2),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    idx_x_scale(calendar, resolved_axis) +
    idx_x_theme(calendar, resolved_axis)
  
  idx_add_info_box(p, calendar, resolved_axis)
}