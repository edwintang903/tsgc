utils::globalVariables(c("Date", "Rt", "lower", "upper", "forecast", "model", "x", "Centered.MA"))

#' @importFrom utils tail
NULL

#' @title Convert an \code{xts}/\code{zoo} object to an \code{idx_series}
#'
#' @description Converts a calendar-indexed \code{xts} or \code{zoo} object
#' with a regular daily frequency (no gaps) into an \code{\link{idx_series}},
#' together with an \code{\link{idx_calendar}} describing how to translate
#' back to calendar time. This is the standard entry point for bringing
#' calendar-indexed data (the user's own \code{xts}/\code{zoo}/data frame
#' data, or one of the package's bundled datasets) into the
#' \code{idx_series}-based analysis functions.
#'
#' @param x An \code{xts} or \code{zoo} object with a daily, gap-free index.
#' @param start.pos Integer position that the first row of \code{x} should
#' be assigned. Use this to align \code{x} with another, already-converted
#' \code{idx_series} that starts at a different calendar date - e.g.
#' \code{start.pos = idx_to_pos(other_cal, zoo::index(x)[1])}. Defaults to
#' \code{1L}.
#'
#' @returns A list with two elements: \code{series}, an \code{idx_series}
#' holding \code{x}'s values, and \code{calendar}, an \code{idx_calendar}
#' anchoring \code{series}'s positions to \code{x}'s original dates.
#'
#' @examples
#' x <- xts::xts(cumsum(rpois(30, 5)) + 1, order.by = Sys.Date() - 29:0)
#' conv <- xts_to_idx(x)
#' conv$series
#' conv$calendar
#'
#' @importFrom xts xts
#' @importFrom zoo index coredata
#'
#' @export
xts_to_idx <- function(x, start.pos = 1L) {
  list(
    series = idx_series(coredata(x), start = start.pos),
    calendar = idx_calendar(
      anchor = index(x)[1],
      anchor_pos = start.pos,
      amount = 1, unit = "days",
      posixct = TRUE
    )
  )
}

#' @title Translate a calendar date to an \code{idx_series} integer position
#'
#' @description Inverse of \code{\link{idx_to_date}}'s calendar-aware
#' branch: converts a calendar date into the integer position (relative to
#' \code{cal}) that corresponds to it. Intended for locating estimation
#' start/end dates, or for aligning a second series to an already-converted
#' \code{idx_series} (see \code{\link{xts_to_idx}}). Requires a
#' calendar-anchored \code{cal} (\code{cal$posixct = TRUE} and a
#' \code{Date}/\code{POSIXct} \code{cal$anchor}); use
#' \code{\link{idx_offset_to_pos}} instead for non-calendar (e.g.
#' sub-second/arbitrary numeric) calendars.
#'
#' For the calendar-relative units \code{"months"}, \code{"quarters"},
#' \code{"years"}, the position is found by counting whole calendar steps
#' between \code{cal$anchor} and \code{date} (the inverse of the
#' calendar-aware stepping used by \code{\link{idx_to_date}}), rather than
#' via a fixed-length difference - a quarter is not a fixed number of days.
#' For the fixed-duration units \code{"seconds"}, \code{"minutes"},
#' \code{"hours"}, \code{"days"}, \code{"weeks"}, positions are computed as
#' a plain time difference scaled by \code{cal$amount}; this is exact for
#' these units. Either way this ignores any \code{pattern} (i.e. exact only
#' for \code{pattern = 1}; only appropriate for simple, evenly spaced
#' calendars).
#'
#' @param cal An \code{idx_calendar} object with \code{posixct = TRUE} and a
#' \code{Date}/\code{POSIXct} \code{anchor}.
#' @param date A single date/time. For \code{cal$unit} of \code{"months"},
#' \code{"quarters"}, or \code{"years"}, coerced via \code{as.Date}. For all
#' other units, coerced via \code{as.POSIXct} if \code{cal$anchor} is
#' \code{POSIXct}, or \code{as.Date} otherwise - so a character string with
#' a time-of-day component (e.g. \code{"2024-01-01 13:30:00"}) is preserved
#' where relevant (sub-day units).
#'
#' @returns A single integer: the \code{idx_series} position corresponding
#' to \code{date}.
#'
#' @examples
#' cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
#'                      amount = 1, unit = "days", posixct = TRUE)
#' idx_to_pos(cal, "2024-01-10")
#'
#' @export
idx_to_pos <- function(cal, date) {
  calendar_step_unit <- cal$unit %in% c("months", "quarters", "years")
  is_calendar_anchor <- inherits(cal$anchor, "Date") || inherits(cal$anchor, "POSIXct")
  
  if (!isTRUE(cal$posixct) || !is_calendar_anchor) {
    stop("idx_to_pos: cal is not a calendar-anchored (posixct) idx_calendar ",
         "(cal$anchor is not a Date/POSIXct, or cal$posixct is FALSE). ",
         "For non-calendar anchors, use idx_offset_to_pos() with a plain ",
         "numeric offset instead of a date.")
  }
  
  if (calendar_step_unit) {
    date <- as.Date(date)
    anchor_date <- as.Date(cal$anchor)
    per_year <- switch(cal$unit, months = 12, quarters = 4, years = 1)
    anchor_ym <- (as.integer(format(anchor_date, "%Y")) * 12 +
                    (as.integer(format(anchor_date, "%m")) - 1))
    date_ym <- (as.integer(format(date, "%Y")) * 12 +
                  (as.integer(format(date, "%m")) - 1))
    months_diff <- date_ym - anchor_ym
    steps <- months_diff / (12 / per_year)
    if (!isTRUE(all.equal(steps, round(steps)))) {
      stop("idx_to_pos: date does not fall on a whole ", cal$unit,
           " boundary relative to cal$anchor.")
    }
    as.integer(cal$anchor_pos + round(steps / cal$amount))
  } else if (inherits(cal$anchor, "POSIXct")) {
    date <- as.POSIXct(date, tz = format(cal$anchor, "%Z"))
    diff_secs <- as.numeric(difftime(date, cal$anchor, units = "secs"))
    by_unit <- idx_calendar_by_unit(cal$unit)
    unit_secs <- switch(if (is.null(by_unit)) "" else by_unit,
                        sec = 1, min = 60, hour = 3600, day = 86400, week = 86400 * 7,
                        NA_real_
    )
    if (is.na(unit_secs)) {
      stop("idx_to_pos: unit '", cal$unit, "' is not a recognised ",
           "fixed-duration unit for a POSIXct anchor.")
    }
    as.integer(cal$anchor_pos + round(diff_secs / unit_secs / cal$amount))
  } else {
    date <- as.Date(date)
    as.integer(cal$anchor_pos + as.numeric(date - cal$anchor) / cal$amount)
  }
}

#' @title Translate a plain numeric offset to an \code{idx_series} integer
#' position
#'
#' @description Inverse of \code{\link{idx_to_date}}'s plain-arithmetic
#' branch (\code{posixct = FALSE}, or a non-Date/POSIXct \code{cal$anchor}):
#' converts a value already expressed in \code{cal$anchor}'s own units/scale
#' (e.g. a picosecond count) into the integer \code{idx_series} position
#' that corresponds to it.
#'
#' Only relevant when you actually have an \code{idx_calendar} with a
#' non-Date/POSIXct \code{anchor} (e.g. a numeric anchor for sub-second or
#' otherwise non-calendar timesteps) and want to look up a position from a
#' value in that anchor's scale. If your series has no calendar meaning at
#' all - no dates, no anchor, nothing to translate - you don't need an
#' \code{idx_calendar} or this function: \code{idx_series} positions
#' (\code{start}/\code{end} arguments throughout this package) are already
#' plain integers, so just use them directly (see \code{\link{idx_range}}).
#'
#' @param cal An \code{idx_calendar} object.
#' @param value A single number, in the same units as \code{cal$anchor}
#' (i.e. \code{cal$anchor + n * cal$amount} for some integer \code{n}), or
#' a vector of such.
#'
#' @returns An integer (vector, matching \code{value}): the \code{idx_series}
#' position(s) corresponding to \code{value}.
#'
#' @examples
#' cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
#'                      unit = "picoseconds")
#' idx_offset_to_pos(cal, 12.5)
#'
#' @export
idx_offset_to_pos <- function(cal, value) {
  stopifnot(is_idx_calendar(cal))
  steps <- (value - cal$anchor) / cal$amount
  if (!isTRUE(all.equal(steps, round(steps)))) {
    stop("idx_offset_to_pos: value does not fall on a whole step boundary ",
         "relative to cal$anchor.")
  }
  as.integer(cal$anchor_pos + round(steps))
}

#' @title Compute log growth rate of cumulated dataset
#
#' @description Helper method to compute the log growth rates of cumulated
#' variables.
#'
#' @param dt Cumulated data series, as an \code{idx_series} with exactly one
#' column.
#' @returns An \code{idx_series} of log growth rates of the cumulated
#' variable inputted via the parameter \code{dt}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1)
#' df2ldl(x)
#'
#' @export
df2ldl <- function(dt) {
  if (!is_idx_series(dt)) {
    stop("dt must be an idx_series object.")
  }
  if (idx_ncol(dt) != 1) {
    stop("dt must only contain 1 data column.")
  }
  lagged <- idx_lag(dt, 1L)
  pos <- idx_positions(dt)
  overlap <- intersect(pos, idx_positions(lagged))
  lag.ov <- lagged[overlap]
  if (any(idx_values(lag.ov) < 0, na.rm = TRUE)) {
    stop("Dataset dt contains negative values.")
  }
  d <- idx_diff(dt, 1L)
  if (any(idx_values(d) < 0, na.rm = TRUE)) {
    stop("Dataset dt has nonpositive increments.")
  }
  lag.aligned <- lagged[idx_positions(d)]
  ldl <- log(idx_values(d) / idx_values(lag.aligned))
  # Pad with a leading NA at dt's start position, since the first
  # observation has no preceding value to compute a growth rate from.
  idx_series(c(NA_real_, ldl), start = dt$start)
}

#' @title Subsetting \code{idx_series} objects given start and end positions
#
#' @description Helper method to subset an \code{idx_series} for a
#' specified range of integer positions.
#'
#' @param df An \code{idx_series} object, or \code{NULL} if no data are
#' supplied.
#' @param start Start position (integer) of the range.
#' @param end End position (integer) of the range. Defaults to the last
#' position in \code{df}.
#' @returns An \code{idx_series} containing the selected observations, or
#'   \code{NULL} if \code{df} is \code{NULL}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1)
#' get_timeframe(x, 5, 10)
#' get_timeframe(x, 5)
#'
#' @export
get_timeframe <- function(df, start, end = NULL) {
  if (is.null(df)) {
    return(NULL)
  }
  if (!is_idx_series(df)) {
    stop("df is not an idx_series object.")
  }
  rng <- idx_range(df)
  if (is.null(end)) {
    end <- rng[2]
  }
  if (length(start) != 1 || length(end) != 1) {
    stop("start and end must each be a single integer position.")
  }
  start <- max(start, rng[1])
  end <- min(end, rng[2])
  if (start > end) {
    stop("start is after end within the available range of df.")
  }
  df[start:end]
}

#' @title Compute successive increments and log growth rate of 2-variable
#' cumulated dataset
#
#' @description Helper method to compute the successive increments and log
#' growth rates of cumulated variables. It will compute the successive
#' increments and log cumulative growth rate for each column in the
#' 2-column series, which will then be used to predict or estimate with the
#' leading indicator model.
#'
#' @param data Cumulated data series as an \code{idx_series} with 2 columns:
#' leading indicator and target variable. Can specify which column is
#' leading indicator by \code{LeadIndCol} parameter.
#' @param LeadIndCol Column number of \code{data} that contains the leading
#' indicator. An integer that can only take values 1 (by default) or 2.
#' @returns A list of \code{idx_series} with the original cumulative
#' variables, successive increments and log growth rates: \code{cLead},
#' \code{cTarg}, \code{newLead}, \code{newTarg}, \code{LDLlead},
#' \code{LDLtarg}.
#'
#' @export
add_daily_ldl <- function(data, LeadIndCol = 1) {
  if (!is_idx_series(data)) {
    stop("data is not an idx_series object.")
  }
  if (idx_ncol(data) != 2) {
    stop("Dataset data must contain exactly two series.")
  }
  mat <- idx_values(data)
  if (LeadIndCol == 1) {
    cLead <- idx_series(mat[, 1], start = data$start)
    cTarg <- idx_series(mat[, 2], start = data$start)
  } else if (LeadIndCol == 2) {
    cLead <- idx_series(mat[, 2], start = data$start)
    cTarg <- idx_series(mat[, 1], start = data$start)
  } else {
    stop("LeadIndCol must be an integer, either 1 or 2.")
  }
  
  newLead <- idx_diff(cLead, 1L)
  newTarg <- idx_diff(cTarg, 1L)
  LDLlead <- df2ldl(cLead)
  LDLtarg <- df2ldl(cTarg)
  
  list(
    cLead = cLead, cTarg = cTarg,
    newLead = newLead, newTarg = newTarg,
    LDLlead = LDLlead, LDLtarg = LDLtarg
  )
}

#' @title Reinitialise a series at a given position
#'
#' @description Takes a cumulated series and re-bases it so that it starts
#' from zero at \code{reinit.idx - 1}.
#'
#' @param dt Cumulated data series, as an \code{idx_series} with exactly one
#' column.
#' @param reinit.idx Integer position at which reinitialisation should
#' occur (i.e. \eqn{t=r}, using the notation in the vignette).
#'
#' @returns The reinitialised series, as an \code{idx_series} starting at
#' \code{reinit.idx}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1)
#' reinitialise_dataframe(x, 10)
#'
#' @export
reinitialise_dataframe <- function(dt, reinit.idx) {
  if (!is_idx_series(dt)) {
    stop("dt is not an idx_series object.")
  }
  if (idx_ncol(dt) != 1) {
    stop("dt must only contain 1 data column.")
  }
  rng <- idx_range(dt)
  if (reinit.idx < rng[1] + 1 || reinit.idx > rng[2]) {
    stop("reinit.idx is not present in dt (or has no preceding value).")
  }
  base <- as.numeric(idx_values(dt[reinit.idx - 1]))
  sub <- dt[reinit.idx:rng[2]]
  idx_series(as.numeric(idx_values(sub)) - base, start = reinit.idx)
}

#' @title Return index and value of maximum
#' @description Similar to Python's argmax function.
#' @param x Object to have its maximum found; either an \code{idx_series}
#' or a plain numeric vector.
#' @param decreasing Logical value indicating whether \code{x} should be
#' ordered in decreasing order. Default is \code{TRUE}. Setting this to
#' \code{FALSE} would find the minimum.
#' @returns If \code{x} is an \code{idx_series}, a length-1 \code{idx_series}
#' at the position of the maximum. Otherwise the maximum value.
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1)
#' argmax(x)
#' @export
argmax <- function(x, decreasing = TRUE) {
  if (is_idx_series(x)) {
    vals <- idx_values(x)
    ord <- order(vals, decreasing = decreasing)[1]
    return(x[idx_positions(x)[ord]])
  }
  return(x[order(x, decreasing = decreasing)[1]])
}

#' @title Compute Mean Absolute Percentage Error (MAPE) for Forecasts Against
#' a Holdout Sample
#'
#' @description This is a helper function that calculates five error metrics
#' of a forecast generated by time series growth curve (tsgc) models. It
#' compares the forecasted values to a holdout sample, providing a measure
#' of forecast accuracy.
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param n.ahead Integer specifying the number of periods to forecast
#' ahead.
#' @param Y An \code{idx_series} object containing the original cumulative
#' dataset.
#'
#' @returns A list containing five error metrics for the forecast, with
#' element names
#' \itemize{
#' \item mape: mean absolute percentage error
#' \item smape: symmetric mean absolute percentage error (between 0 to 100)
#' \item mae: mean absolute error
#' \item rmse: root mean squared error
#' \item coverage: Percentage of holdout sample data points that lie inside
#' the confidence interval for predictions}
#'
#' @export
mapes <- function(res, n.ahead, Y) {
  res$mapes(n.ahead, Y)
}

#' @title Estimate the reproduction number from a fitted model
#'
#' @description Computes the (instantaneous) reproduction number
#' \eqn{R_t = \exp\{g_t \times \code{gen_int}\}}, where \eqn{g_t} is the
#' filtered or smoothed growth rate of the incidence variable \eqn{y}
#' returned by \code{res$get_gy_ci()}, and \code{gen_int} is the mean
#' generation interval of the disease/process being modelled (in the same
#' integer-position units as \code{res}). This estimate, and its confidence
#' interval, is a deterministic transformation of \code{get_gy_ci()}'s
#' output; \code{res} must already have been produced by \code{estimate()}.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param gen_int The mean generation interval, in the same integer-position
#' units as \code{res} (e.g. if positions are days, this is the mean
#' generation interval in days).
#' @param n.ahead Number of most recent integer positions to return, taken
#' from the end of \code{res}'s filtered/smoothed range. Default is
#' \code{7}.
#' @param smoothed Logical value indicating whether to use the smoothed
#' (\code{TRUE}) or filtered (\code{FALSE}, default) growth rate estimates
#' underlying \eqn{R_t}. Passed through to \code{res$get_gy_ci()}.
#' @param confidence.level Confidence level for the confidence interval
#' around \eqn{R_t}. Default is \code{0.68}, one standard deviation for a
#' normally distributed random variable. Passed through to
#' \code{res$get_gy_ci()}.
#'
#' @returns An \code{idx_series} with columns \code{fit}, \code{lower} and
#' \code{upper}, giving the estimated reproduction number and its
#' confidence interval over the last \code{n.ahead} integer positions.
#'
#' @examples
#' library(tsgc)
#' set.seed(1)
#' Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1)
#' model <- SSModelDynamicGompertz$new(Y = Y, q = NULL, end = 100)
#' res <- estimate(model)
#' estimate_r0(res, gen_int = 5, n.ahead = 7)
#'
#' @export
estimate_r0 <- function(res, gen_int, n.ahead = 7, smoothed = FALSE,
                        confidence.level = 0.68) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  gy.ci <- res$get_gy_ci(smoothed = smoothed, confidence.level = confidence.level)
  rt_mat <- exp(idx_values(gy.ci) * gen_int)
  n <- nrow(rt_mat)
  keep <- seq.int(max(1L, n - n.ahead + 1L), n)
  idx_series(rt_mat[keep, , drop = FALSE], start = gy.ci$start + keep[1] - 1L)
}

#' @title Walk-Forward Validation for Model Comparison Using Mean Absolute
#' Percentage Error (MAPE)
#'
#' @description This function performs a walk-forward validation to compare
#' forecasting performance across different models specified by the user.
#' It returns a data frame of a user-specified error metric (e.g. MAPE, MAE)
#' for forecasts \code{n.ahead} positions ahead, using the given models with
#' varying end positions.
#'
#' @param Y An \code{idx_series} representing the cumulative data series. If
#' a Leading Indicator model is compared, Y should include columns for both
#' the leading indicator and the target variable. The specific column for
#' the leading indicator can be designated using the \code{LeadIndCol}
#' parameter.
#' @param model_list A list containing \code{SSModelDynamicGompertz} or
#' \code{SSModelLeadingIndicator} objects, to be compared in a cross
#' validation procedure.
#' @param est.end The initial estimation end position for model fitting.
#' Starting from this position, the function re-estimates the model and
#' evaluates the performance for each lag in \code{all_lags} every
#' \code{gap} positions, over \code{n.estimate} steps.
#' @param n.ahead Integer specifying the number of positions to forecast
#' ahead for MAPE evaluation.
#' @param n.estimate Integer indicating the total number of walk-forward
#' validation steps to report.
#' @param gap Integer specifying the position gap between two successive
#' validations, where the model is re-estimated and evaluated during the
#' walk-forward validation.
#' @param xpred_lead.full (Only required for leading indicator models) An
#' \code{idx_series} containing the values of exogenous variables for the
#' leading indicator over the estimation and prediction time frame.
#' @param xpred_targ.full An \code{idx_series} containing the values of
#' exogenous variables for the target variable over the estimation and
#' prediction time frame.
#' @param LeadIndCol (Only required for leading indicator models) Integer
#' representing the column number in \code{Y} that contains the leading
#' indicator.
#' @param criterion A string object indicating how to compare between
#' different models. Available choices are "mape" (by default), "smape",
#' "mae" and "rmse".
#'
#' @returns A table summarizing the chosen error metric for each model in
#' \code{model_list} across the specified positions.
#'
#' @export
cross_val <- function(Y, model_list, est.end, n.ahead = 7, n.estimate = 1, gap = 1,
                      xpred_targ.full = NULL, xpred_lead.full = NULL,
                      LeadIndCol = 1, criterion = "mape") {
  if (!is_idx_series(Y)) {
    stop("Y must be an idx_series object.")
  }
  if (idx_ncol(Y) == 1) {
    Y1 <- Y
  } else if (idx_ncol(Y) == 2) {
    keep_col <- setdiff(1:2, LeadIndCol)
    Y1 <- idx_series(idx_values(Y)[, keep_col], start = Y$start)
  } else {
    stop("Y should not have more than 2 columns.")
  }
  if (length(est.end) != 1 || !isTRUE(all.equal(est.end, as.integer(est.end)))) {
    stop("est.end must be a single integer position.")
  }
  if (n.ahead <= 0) {
    stop("n.ahead must be a positive integer.")
  }
  results <- data.frame(
    Model = names(model_list)
  )
  for (k in 1:n.estimate) {
    index_num <- 1
    for (model in model_list) {
      model$end <- est.end + (k - 1) * gap
      if (inherits(model, "SSModelDynamicGompertz")) {
        model$Y <- get_timeframe(Y1, model$start, model$end)
        if (!is.null(model$xpred)) {
          model$xpred <- get_timeframe(xpred_targ.full, model$start, model$end)
        }
        res <- estimate(model)
        if (res$xpred_logical) {
          res$xpred.new <- xpred_targ.full
        }
        results[index_num, k + 1] <- round(mapes(res, n.ahead, Y1)[[criterion]], 2)
      } else if (inherits(model, "SSModelLeadingIndicator")) {
        if (!is.null(model$xpred_lead)) {
          model$xpred_lead <- xpred_lead.full
        }
        if (!is.null(model$xpred_targ)) {
          model$xpred_targ <- xpred_targ.full
        }
        res <- estimate(model)
        if (res$xpred_logical[1]) {
          res$xpred_lead.new <- xpred_lead.full
        }
        if (res$xpred_logical[2]) {
          res$xpred_targ.new <- xpred_targ.full
        }
        results[index_num, k + 1] <- round(mapes(res, n.ahead, Y)[[criterion]], 2)
      } else {
        stop(paste("Model", index_num, "in model_list is not a SSModelDynamicGompertz or SSModelLeadingIndicator object."))
      }
      index_num <- index_num + 1
    }
  }
  all_ends <- as.character(est.end + c(0:(k - 1)) * gap)
  colnames(results) <- c("Model", all_ends)
  return(results)
}

#' @title Write a selection of relevant results to disk
#'
#' @description Function writes the following results to csv files which get
#' saved in the location specified in \code{res.dir}: forecast new cases or
#' incidence variable, \eqn{y}; the filtered level and slope of \eqn{\ln g},
#' \eqn{\delta} and \eqn{\gamma}; filtered estimates of \eqn{g_y} and the
#' confidence intervals for these estimates.
#'
#' @param res Results object of class \code{FilterResults} or \code{FilterResultsLI},
#' obtained from \samp{estimate()} method.
#' @param res.dir File path to save the results to. A character string.
#' @param n.ahead Number of periods ahead to forecast. A positive integer.
#' @param prefix The prefix to be added to the file names generated. A character string.
#' @param confidence.level Confidence level to use for the confidence interval
#' on the forecasts \eqn{\ln(g_t)}.
#'
#' @importFrom utils write.csv
#' @importFrom stats qnorm
#'
#' @returns A number of csv files saved in the directory specified in
#' \code{res.dir}.
#' @examples
#' # Not run as do not wish to save to local disk when compiling documentation.
#' # Below will run if copied and pasted into console.
#' library(tsgc)
#' library(here)
#'
#' res.dir <- tempdir()
#' data(gauteng,package="tsgc")
#' conv <- xts_to_idx(gauteng)
#' res <- estimate(SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#' calendar = conv$calendar,
#' end = idx_to_pos(conv$calendar, "2020-07-06")))
#'
#' tsgc::write_results(
#' res=res, res.dir = res.dir, prefix="dyn_gompertz",n.ahead = 14,
#' confidence.level = 0.68)
#'
#' @export
write_results <- function(res, res.dir, n.ahead, prefix="", confidence.level=0.68) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")){
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  
  # New cases (delta Y)
  y.hat.diff <- res$predict_level(
    n.ahead = n.ahead,
    confidence.level = confidence.level,
    sea.on = TRUE)
  
  write.csv(
    as.matrix(y.hat.diff),
    row.names = idx_positions(y.hat.diff),
    file = file.path(res.dir, paste(prefix, "cases_fcst.csv", sep="")))
  
  # Filtered slope / level
  y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
  filtered.level <- y.hat.all$level.t.t
  filtered.slope <- y.hat.all$slope.t.t
  a.t.t <- y.hat.all$a.t.t
  P.t.t <- y.hat.all$P.t.t
  idx.slope <- grep("slope", colnames(a.t.t))
  idx.level <- grep("level", colnames(a.t.t))[1]
  gamma.std.err <- sqrt(P.t.t[idx.slope, idx.slope,])
  delta.std.err <- sqrt(P.t.t[idx.level, idx.level,])
  
  gamma <- idx_cbind(filtered.slope, idx_series(gamma.std.err, start = filtered.slope$start))
  delta <- idx_cbind(filtered.level, idx_series(delta.std.err, start = filtered.level$start))
  gamma.mat <- as.matrix(gamma)
  delta.mat <- as.matrix(delta)
  colnames(gamma.mat) <- c("gamma", "std.err")
  colnames(delta.mat) <- c("delta", "std.err")
  
  write.csv(
    gamma.mat,
    row.names = idx_positions(filtered.slope),
    file = file.path(res.dir, paste(prefix, "trend_slope_filt.csv", sep=""))
  )
  write.csv(
    delta.mat,
    row.names = idx_positions(filtered.level),
    file = file.path(res.dir, paste(prefix, "log_gr_level_filt.csv", sep=""))
  )
  
  # Filtered growth rate of new cases (g_y); CI from standard error on the
  # slope component of the state covariance matrix.
  g.y.t.t <- exp(idx_values(filtered.level)) + idx_values(filtered.slope)
  ci <- qnorm((1 - confidence.level) / 2) * gamma.std.err %o% c(1, -1)
  ci_bounds <- as.vector(g.y.t.t) + ci
  gy.ci.mat <- cbind(fit = g.y.t.t, prediction = ci_bounds)
  colnames(gy.ci.mat) <- c("fit", "lower", "upper")
  
  write.csv(
    gy.ci.mat,
    row.names = idx_positions(filtered.level),
    file = file.path(res.dir, paste(prefix, "cases_gr.csv", sep="")))
  
  message("Saved results for: ", substitute(res))
}