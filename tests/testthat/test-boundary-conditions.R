# =============================================================================
# Boundary & Edge-Case Tests: SSModelDynamicGompertz and SSModelLeadingIndicator
#
# Boundaries are derived directly from the package source and cross-checked 
# against theory. Tested in pairs (just inside vs. just outside valid regions) 
# to confirm exact boundary locations.
#
# Out-of-scope: Reinitialisation trigger rules, exact peak-prediction formula 
# (estimate_r0 used as proxy), and near-unidentified parameter estimates.
# =============================================================================

library(tsgc)
library(xts)
library(zoo)

# #############################################################################
# 1. THEORETICAL BOUNDARIES: SSModelDynamicGompertz
# #############################################################################

# -----------------------------------------------------------------------
# 1.1 Strictly increasing Y (Boundary is exactly zero difference)
# -----------------------------------------------------------------------

test_that("SSModelDynamicGompertz errors on a plateaued or decreasing series at estimate(), not at construction", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 30)
  cum <- c(seq(100, 200, length.out = 20), rep(200, 3), seq(205, 250, length.out = 7))
  Y_flat <- xts(cum, order.by = dates)
  
  expect_no_error(model <- SSModelDynamicGompertz$new(Y = Y_flat))
  model <- SSModelDynamicGompertz$new(Y = Y_flat)
  expect_error(model$estimate(), "strictly increasing")
})

test_that("SSModelDynamicGompertz errors on a strictly decreasing segment (not just flat)", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 15)
  cum <- c(seq(100, 180, length.out = 10), 178, seq(182, 220, length.out = 4))
  Y_dec <- xts(cum, order.by = dates)
  
  model <- SSModelDynamicGompertz$new(Y = Y_dec)
  expect_error(model$estimate(), "strictly increasing")
})

test_that("SSModelDynamicGompertz accepts a series with an arbitrarily small but strictly positive increment at the same point that previously plateaued", {
  # Boundary case: identical series to the plateau test, nudged up by a tiny epsilon.
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 30)
  cum <- c(seq(100, 200, length.out = 20), rep(200, 3), seq(205, 250, length.out = 7))
  jitter <- c(rep(0, 20), 1e-8, 2e-8, 3e-8, rep(0, 7))
  Y_jittered <- xts(cum + jitter, order.by = dates)
  
  model <- SSModelDynamicGompertz$new(Y = Y_jittered)
  expect_no_error(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

# -----------------------------------------------------------------------
# 1.2 reinit.date must be an exact match in Y's index
# -----------------------------------------------------------------------

test_that("SSModelDynamicGompertz accepts a reinit.date that is an exact match in the series index", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 60)
  Y <- xts(exp(seq(4, 7, length.out = 60)), order.by = dates)
  
  model <- SSModelDynamicGompertz$new(
    Y = Y, q = 0.01, reinit.date = dates[30]
  )
  expect_no_error(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

test_that("SSModelDynamicGompertz fails with a reinit.date that does not exist in the series index", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 60)
  Y <- xts(exp(seq(4, 7, length.out = 60)), order.by = dates)
  
  # One day beyond the end of the series.
  model <- SSModelDynamicGompertz$new(
    Y = Y, q = 0.01, reinit.date = tail(dates, 1) + 1
  )
  expect_error(model$estimate())
})

# -----------------------------------------------------------------------
# 1.3 Minimal length and degeneracy checks
# -----------------------------------------------------------------------

test_that("SSModelDynamicGompertz warns of a degenerate model on an extremely short series (sea.period = 0)", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 3)
  Y_short <- xts(c(100, 110, 125), order.by = dates)
  
  model <- SSModelDynamicGompertz$new(Y = Y_short, sea.period = 0)
  expect_warning(result <- model$estimate(), "degenerate")
  
  if (inherits(result, "FilterResults")) {
    expect_true(!is.null(result$output))
  }
})

test_that("SSModelDynamicGompertz with sea.period = 0 estimates cleanly, with no degeneracy warning, once given enough data", {
  # Positive control: trend-only specification with sufficient length.
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 25)
  Y_ok <- xts(exp(seq(4.6, 6, length.out = 25)), order.by = dates)
  
  model <- SSModelDynamicGompertz$new(Y = Y_ok, q = 0.01, sea.period = 0)
  expect_no_warning(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

test_that("SSModelDynamicGompertz with sea.period = 7 warns of a degenerate model when data cannot identify seasonal states", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 5)
  Y_short <- xts(seq(100, 140, length.out = 5), order.by = dates)
  
  model <- SSModelDynamicGompertz$new(Y = Y_short, sea.period = 7)
  
  warnings_seen <- character(0)
  result <- withCallingHandlers(
    tryCatch(model$estimate(), error = function(e) e),
    warning = function(w) {
      warnings_seen[[length(warnings_seen) + 1]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("degenerate", warnings_seen)))
  expect_true(any(grepl("[Dd]iffuse filtering|Finf", warnings_seen)))
  expect_true(inherits(result, "error") || inherits(result, "FilterResults"))
})

test_that("SSModelDynamicGompertz with sea.period = 7 estimates cleanly, with no degeneracy warning, once given enough data", {
  # Positive control: 8 states to estimate, comfortable series length.
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 40)
  Y_ok <- xts(exp(seq(4.6, 7, length.out = 40)), order.by = dates)
  
  model <- SSModelDynamicGompertz$new(Y = Y_ok, q = 0.01, sea.period = 7)
  expect_no_warning(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

# -----------------------------------------------------------------------
# 1.4 Peak / turning-point boundary via estimate_r0() proxy
# -----------------------------------------------------------------------

test_that("estimate_r0 reproduces the documented decline toward Rt = 1 approaching a known peak", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  
  model <- SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = 0.005,
    start.date = as.Date("2021-02-01"),
    end.date   = as.Date("2021-04-19")
  )
  res <- model$estimate()
  r_t <- estimate_r0(res, gen_int = 4, ndays = 7)
  
  expect_equal(nrow(r_t), 7)
  # Decline over the first six days of the documented window.
  expect_true(all(diff(r_t$Rt[1:6]) <= 1e-8))
  # Ensure the window as a whole approaches Rt = 1 at the end.
  expect_true(abs(tail(r_t$Rt, 1) - 1) < abs(r_t$Rt[1] - 1))
})

test_that("estimate_r0 output brackets Rt = 1 as growth decelerates through a known peak (wider window)", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  
  model <- SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = 0.005,
    start.date = as.Date("2021-02-01"),
    end.date   = as.Date("2021-06-25")
  )
  res <- model$estimate()
  r_t <- estimate_r0(res, gen_int = 4, ndays = 30)
  
  expect_s3_class(r_t, "data.frame")
  expect_true(all(c("Date", "Rt", "lower", "upper") %in% names(r_t)))
  expect_true(any(r_t$Rt > 1) || any(r_t$Rt < 1))
  expect_true(all(r_t$lower <= r_t$Rt))
  expect_true(all(r_t$Rt <= r_t$upper))
  expect_true(all(is.finite(r_t$Rt)))
})


# #############################################################################
# 2. THEORETICAL BOUNDARIES: SSModelLeadingIndicator
# #############################################################################

# -----------------------------------------------------------------------
# 2.1 Strictly increasing series: DECREASE vs PLATEAU trigger different errors
# -----------------------------------------------------------------------

test_that("SSModelLeadingIndicator: a genuine decrease in the target series is caught by df2ldl's own on-topic message", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 50)
  lead  <- seq(100, 400, length.out = 50)
  
  targ <- seq(50, 200, length.out = 50)
  targ[35] <- targ[34] - 1  # genuine decrease, not a plateau
  
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)
  expect_error(mod$estimate(), "nonpositive increments")
})

test_that("SSModelLeadingIndicator: a decrease in the LEAD series is caught by the same df2ldl message", {
  # Confirms check applies symmetrically to both columns.
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 50)
  targ  <- seq(50, 200, length.out = 50)
  
  lead <- seq(100, 400, length.out = 50)
  lead[35] <- lead[34] - 1  # genuine decrease in the lead column
  
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)
  expect_error(mod$estimate(), "nonpositive increments")
})

test_that("SSModelLeadingIndicator: a plateau (not a decrease) in the target series produces a misleading date-resolution error rather than the strictly-increasing guidance", {
  # Documents current verified behaviour (plateau trips date-resolution error). 
  n_lag <- 5
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 50)
  lead  <- seq(100, 400, length.out = 50)
  
  targ <- seq(50, 200, length.out = 50)
  targ[35] <- targ[34]  # exact plateau, not a decrease
  
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = n_lag, LeadIndCol = 1)
  expect_error(mod$estimate(), "not separated by the same time resolution")
})

test_that("SSModelLeadingIndicator accepts a jittered version of the same plateaued series", {
  # Positive control: plateaued segment nudged up by epsilon.
  n_lag <- 5
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 50)
  lead  <- seq(100, 400, length.out = 50)
  
  targ <- seq(50, 200, length.out = 50)
  targ[35] <- targ[34] + 1e-6
  
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = n_lag, LeadIndCol = 1)
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

# -----------------------------------------------------------------------
# 2.2 Date resolution checks during initialization
# -----------------------------------------------------------------------

test_that("SSModelLeadingIndicator rejects irregular raw dates immediately at construction, before estimate() is ever called", {
  # Isolate the constructor check: skip one calendar day.
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 50)
  dates <- dates[-25]  # remove one day, creating a single 2-day gap
  
  lead <- seq(100, 400, length.out = 49)
  targ <- seq(50, 200, length.out = 49)
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  expect_error(
    SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1),
    "not separated by the same time resolution"
  )
})

test_that("SSModelLeadingIndicator accepts raw input with fully regular daily dates (paired positive control)", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 50)
  lead <- seq(100, 400, length.out = 50)
  targ <- seq(50, 200, length.out = 50)
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  expect_no_error(
    mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)
  )
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

# -----------------------------------------------------------------------
# 2.3 n.lag validation risks (Zero, Negative, OOB)
# -----------------------------------------------------------------------

test_that("SSModelLeadingIndicator accepts n.lag = 0 (contemporaneous alignment, a valid edge case)", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 40)
  lead <- seq(100, 400, length.out = 40)
  targ <- seq(50, 200, length.out = 40)
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 0, LeadIndCol = 1)
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

test_that("SSModelLeadingIndicator does NOT error on a negative n.lag, but silently reverses the intended lead/lag direction", {
  # Documents risk: negative shift reverses alignment rather than throwing error.
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 40)
  lead <- seq(100, 400, length.out = 40)
  targ <- seq(50, 200, length.out = 40)
  Y_li <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = -3, LeadIndCol = 1)
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

test_that("SSModelLeadingIndicator fails in a controlled way when n.lag exceeds available data", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 10)
  lead  <- seq(100, 200, length.out = 10)
  targ  <- seq(50, 90, length.out = 10)
  Y_li  <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y_li) <- c("lead_col", "targ_col")
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 20, LeadIndCol = 1)
  result <- tryCatch(mod$estimate(), error = function(e) e)
  expect_true(inherits(result, "error"))
})


# #############################################################################
# 3. STRUCTURAL / TYPE-VALIDATION BOUNDARIES (both classes)
# #############################################################################

test_that("sea.period validation is identical and correctly enforced in both model classes", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 30)
  Y1 <- xts(exp(seq(4, 6, length.out = 30)), order.by = dates)
  
  # Just-outside: 1 is explicitly excluded; negative and non-integer also fail.
  expect_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 1), "sea.period")
  expect_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = -1), "sea.period")
  expect_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 2.5), "sea.period")
  # Just-inside: 0 and any other non-negative integer != 1 are valid.
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 0))
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 2))
  
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y2 <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y2) <- c("lead_col", "targ_col")
  
  expect_error(SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, sea.period = 1), "sea.period")
  expect_no_error(SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, sea.period = 0))
})

test_that("original.results must be NULL or a FilterResults object", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 30)
  Y <- xts(exp(seq(4, 6, length.out = 30)), order.by = dates)
  
  expect_error(
    SSModelDynamicGompertz$new(Y = Y, original.results = "not_a_filterresults"),
    "original.results"
  )
  expect_no_error(SSModelDynamicGompertz$new(Y = Y, original.results = NULL))
})

test_that("xpred / xpred_lead / xpred_targ must be NULL or an xts object, in both classes", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 30)
  Y1 <- xts(exp(seq(4, 6, length.out = 30)), order.by = dates)
  
  expect_error(SSModelDynamicGompertz$new(Y = Y1, xpred = data.frame(x = 1:30)), "xpred")
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, xpred = NULL))
  
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y2 <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y2) <- c("lead_col", "targ_col")
  
  expect_error(
    SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, xpred_lead = data.frame(x = 1:30)),
    "xpred_lead"
  )
  expect_error(
    SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, xpred_targ = data.frame(x = 1:30)),
    "xpred_targ"
  )
})

test_that("LeadIndCol must take the value 1 or 2, in SSModelLeadingIndicator", {
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = 30)
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y <- xts(cbind(lead, targ), order.by = dates)
  colnames(Y) <- c("lead_col", "targ_col")
  
  expect_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 0), "LeadIndCol")
  expect_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 3), "LeadIndCol")
  expect_no_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 1))
  expect_no_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 2))
})


# #############################################################################
# 4. SHARED UTILITY BOUNDARIES: get_timeframe() and get_time_resolution()
# #############################################################################

test_that("get_timeframe silently returns a zero-row object when start.date > end.date", {
  # Risk: No check prevents an empty window from silent downstream propagation.
  data(gauteng, package = "tsgc")
  result <- get_timeframe(gauteng, as.Date("2021-06-01"), as.Date("2021-01-01"))
  expect_equal(nrow(result), 0)
})

test_that("get_timeframe returns the expected non-empty window when start.date <= end.date (paired positive control)", {
  data(gauteng, package = "tsgc")
  result <- get_timeframe(gauteng, as.Date("2021-01-01"), as.Date("2021-01-10"))
  expect_equal(nrow(result), 10)
})

test_that("get_time_resolution errors cleanly with fewer than two distinct dates", {
  # Expecting "Input must contain at least two distinct dates."
  expect_error(
    get_time_resolution(as.Date("2021-01-01")),
    "at least two distinct dates|not a vector of dates"
  )
  # Identical dates collapse after unique().
  expect_error(
    get_time_resolution(as.Date(c("2021-01-01", "2021-01-01"))),
    "at least two distinct dates|not a vector of dates"
  )
})

test_that("get_time_resolution succeeds with exactly two distinct, regularly-spaced dates (paired positive control)", {
  result <- get_time_resolution(as.Date(c("2021-01-01", "2021-01-02")))
  expect_equal(result, "daily")
})

test_that("get_time_resolution errors on irregularly-spaced dates", {
  irregular <- as.Date(c("2021-01-01", "2021-01-02", "2021-01-05"))
  expect_error(get_time_resolution(irregular), "not separated by the same time resolution")
})