# Created by: Craig Thamotheram
# Created on: 19/02/2022
# Refactored: analysis functions operate on idx_series (integer-indexed),
# not on calendar time.

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

utils::globalVariables(c("Date", "Rt", "lower", "upper", "forecast", "model"))

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
  idx_series(log(idx_values(d) / idx_values(lag.aligned)), start = d$start)
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
  base <- idx_values(dt[reinit.idx - 1])
  sub <- dt[reinit.idx:rng[2]]
  idx_series(idx_values(sub) - base, start = reinit.idx)
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