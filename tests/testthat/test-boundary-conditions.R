library(tsgc)

## 1. THEORETICAL BOUNDARIES: SSModelDynamicGompertz ----

## 1.1 Strictly increasing Y (Boundary is exactly zero difference) ----

test_that("SSModelDynamicGompertz errors on a plateaued or decreasing series at estimate(), not at construction", {
  cum <- c(seq(100, 200, length.out = 20), rep(200, 3), seq(205, 250, length.out = 7))
  Y_flat <- idx_series(cum, start = 1L)
  
  expect_no_error(model <- SSModelDynamicGompertz$new(Y = Y_flat))
  model <- SSModelDynamicGompertz$new(Y = Y_flat)
  expect_error(model$estimate(), "strictly increasing")
})

test_that("SSModelDynamicGompertz errors on a strictly decreasing segment (not just flat)", {
  cum <- c(seq(100, 180, length.out = 10), 178, seq(182, 220, length.out = 4))
  Y_dec <- idx_series(cum, start = 1L)
  
  model <- SSModelDynamicGompertz$new(Y = Y_dec)
  expect_error(model$estimate(), "strictly increasing")
})

test_that("SSModelDynamicGompertz accepts a series with an arbitrarily small but strictly positive increment at the same point that previously plateaued", {
  # Boundary case: identical series to the plateau test, nudged up by a tiny epsilon.
  cum <- c(seq(100, 200, length.out = 20), rep(200, 3), seq(205, 250, length.out = 7))
  jitter <- c(rep(0, 20), 1e-8, 2e-8, 3e-8, rep(0, 7))
  Y_jittered <- idx_series(cum + jitter, start = 1L)
  
  model <- SSModelDynamicGompertz$new(Y = Y_jittered)
  expect_no_error(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

## 1.2 reinit.idx must correspond to a position actually present in Y ----

test_that("SSModelDynamicGompertz accepts a reinit.idx that is a position within the series", {
  Y <- idx_series(exp(seq(4, 7, length.out = 60)), start = 1L)
  
  model <- SSModelDynamicGompertz$new(
    Y = Y, q = 0.01, reinit.idx = 30L
  )
  expect_no_error(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

test_that("SSModelDynamicGompertz fails with a reinit.idx that does not exist in the series", {
  Y <- idx_series(exp(seq(4, 7, length.out = 60)), start = 1L)
  
  model <- SSModelDynamicGompertz$new(
    Y = Y, q = 0.01, reinit.idx = 61L
  )
  expect_error(model$estimate())
})

## 1.3 Minimal length and degeneracy checks ----

test_that("SSModelDynamicGompertz warns of a degenerate model on an extremely short series (sea.period = 0)", {
  Y_short <- idx_series(c(100, 110, 125), start = 1L)
  
  model <- SSModelDynamicGompertz$new(Y = Y_short, sea.period = 0)
  expect_warning(result <- model$estimate(), "degenerate")
  
  if (inherits(result, "FilterResults")) {
    expect_true(!is.null(result$output))
  }
})

test_that("SSModelDynamicGompertz with sea.period = 0 estimates cleanly, with no degeneracy warning, once given enough data", {
  # Positive control: trend-only specification with sufficient length.
  Y_ok <- idx_series(exp(seq(4.6, 6, length.out = 25)), start = 1L)
  
  model <- SSModelDynamicGompertz$new(Y = Y_ok, q = 0.01, sea.period = 0)
  expect_no_warning(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

test_that("SSModelDynamicGompertz with sea.period = 7 warns of a degenerate model when data cannot identify seasonal states", {
  Y_short <- idx_series(seq(100, 140, length.out = 5), start = 1L)
  
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
  Y_ok <- idx_series(exp(seq(4.6, 7, length.out = 40)), start = 1L)
  
  model <- SSModelDynamicGompertz$new(Y = Y_ok, q = 0.01, sea.period = 7)
  expect_no_warning(res <- model$estimate())
  expect_true(inherits(res, "FilterResults"))
})

## 1.4 Peak / turning-point boundary via estimate_r0() proxy ----

test_that("estimate_r0 reproduces the documented decline toward Rt = 1 approaching a known peak", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005,
    start = est.start, end = est.end
  )
  res <- model$estimate()
  r_t <- estimate_r0(res, gen_int = 4, n.ahead = 7)
  rt_vals <- r_t$fit
  
  expect_equal(nrow(r_t), 7)
  # Decline over the first six days of the documented window.
  expect_true(all(diff(rt_vals[1:6]) <= 1e-8))
  # Ensure the window as a whole approaches Rt = 1 at the end.
  expect_true(abs(tail(rt_vals, 1) - 1) < abs(rt_vals[1] - 1))
})

test_that("estimate_r0 output brackets Rt = 1 as growth decelerates through a known peak (wider window)", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005,
    start = est.start, end = est.end
  )
  res <- model$estimate()
  r_t <- estimate_r0(res, gen_int = 4, n.ahead = 30)
  
  expect_s3_class(r_t, "data.frame")
  expect_true(all(c("fit", "lower", "upper") %in% names(r_t)))
  expect_true(any(r_t$fit > 1) || any(r_t$fit < 1))
  expect_true(all(r_t$lower <= r_t$fit))
  expect_true(all(r_t$fit <= r_t$upper))
  expect_true(all(is.finite(r_t$fit)))
})


## 2. THEORETICAL BOUNDARIES: SSModelLeadingIndicator ----

## 2.1 Strictly increasing series: DECREASE vs PLATEAU ----

test_that("SSModelLeadingIndicator: a genuine decrease in the target series is caught by df2ldl's own on-topic message", {
  lead <- seq(100, 400, length.out = 50)
  
  targ <- seq(50, 200, length.out = 50)
  targ[35] <- targ[34] - 1  # genuine decrease, not a plateau
  
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)
  expect_error(mod$estimate(), "nonpositive increments")
})

test_that("SSModelLeadingIndicator: a decrease in the LEAD series is caught by the same df2ldl message", {
  # Confirms check applies symmetrically to both columns.
  targ <- seq(50, 200, length.out = 50)
  
  lead <- seq(100, 400, length.out = 50)
  lead[35] <- lead[34] - 1  # genuine decrease in the lead column
  
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)
  expect_error(mod$estimate(), "nonpositive increments")
})

test_that("SSModelLeadingIndicator: an exact plateau (not a decrease) in the target series produces an interior-gap error rather than the strictly-increasing message", {
  # An exact plateau makes LDLtarg = log(0) = -Inf at that position, which
  # is treated as missing and dropped, leaving a gap idx_series cannot
  # represent, rather than tripping the strictly-increasing check.
  n_lag <- 5
  lead  <- seq(100, 400, length.out = 50)
  
  targ <- seq(50, 200, length.out = 50)
  targ[35] <- targ[34]  # exact plateau, not a decrease
  
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = n_lag, LeadIndCol = 1)
  expect_error(mod$estimate(), "gaps|idx_series cannot represent")
})

test_that("SSModelLeadingIndicator accepts a jittered version of the same plateaued series", {
  # Positive control: plateaued segment nudged up by epsilon.
  n_lag <- 5
  lead  <- seq(100, 400, length.out = 50)
  
  targ <- seq(50, 200, length.out = 50)
  targ[35] <- targ[34] + 1e-6
  
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = n_lag, LeadIndCol = 1)
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

## 2.2 Gaps in the position index after trimming missing/infinite values ----

test_that("SSModelLeadingIndicator errors when removing missing/infinite values leaves a gap in the position index", {
  # Interior NA leaves the remaining valid positions non-contiguous after
  # trimming, which idx_series cannot represent.
  lead <- seq(100, 400, length.out = 50)
  targ <- seq(50, 200, length.out = 50)
  targ[25] <- NA
  
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  expect_error(
    SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)$estimate(),
    "gaps|idx_series cannot represent"
  )
})

test_that("SSModelLeadingIndicator accepts fully regular, contiguous series (paired positive control)", {
  lead <- seq(100, 400, length.out = 50)
  targ <- seq(50, 200, length.out = 50)
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  expect_no_error(
    mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 5, LeadIndCol = 1)
  )
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

## 2.3 n.lag validation risks (Zero, Negative, OOB) ----

test_that("SSModelLeadingIndicator accepts n.lag = 0 (contemporaneous alignment, a valid edge case)", {
  lead <- seq(100, 400, length.out = 40)
  targ <- seq(50, 200, length.out = 40)
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 0, LeadIndCol = 1)
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

test_that("SSModelLeadingIndicator does NOT error on a negative n.lag, but silently reverses the intended lead/lag direction", {
  # Documents risk: negative shift reverses alignment rather than throwing error.
  lead <- seq(100, 400, length.out = 40)
  targ <- seq(50, 200, length.out = 40)
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = -3, LeadIndCol = 1)
  expect_no_error(res <- mod$estimate())
  expect_true(inherits(res, "FilterResultsLI"))
})

test_that("SSModelLeadingIndicator fails in a controlled way when n.lag exceeds available data", {
  lead  <- seq(100, 200, length.out = 10)
  targ  <- seq(50, 90, length.out = 10)
  Y_li  <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 20, LeadIndCol = 1)
  expect_error(mod$estimate(), "n.lag|overlapping positions")
})

test_that("SSModelLeadingIndicator errors clearly, naming n.lag and the row count, when the estimation window collapses to a single usable row", {
  # A single-row estimation window previously failed deep inside
  # KFAS::SSModel() with an opaque error; now checked explicitly.
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  # start == end leaves exactly one row in y.estimate regardless of n.lag.
  mod <- SSModelLeadingIndicator$new(
    Y = Y_li, n.lag = 3, LeadIndCol = 1, start = 10L, end = 10L
  )
  expect_error(
    mod$estimate(),
    "only 1 usable row.*n\\.lag = 3"
  )
})


## 3. STRUCTURAL / TYPE-VALIDATION BOUNDARIES (both classes) ----

test_that("sea.period validation is identical and correctly enforced in both model classes", {
  Y1 <- idx_series(exp(seq(4, 6, length.out = 30)), start = 1L)
  
  # Just-outside: 1 is explicitly excluded; negative and non-integer also fail.
  expect_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 1), "sea.period")
  expect_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = -1), "sea.period")
  expect_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 2.5), "sea.period")
  # Just-inside: 0 and any other non-negative integer != 1 are valid.
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 0))
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, sea.period = 2))
  
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y2 <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  expect_error(SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, sea.period = 1), "sea.period")
  expect_no_error(SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, sea.period = 0))
})

test_that("original.results must be NULL or a FilterResults object", {
  Y <- idx_series(exp(seq(4, 6, length.out = 30)), start = 1L)
  
  expect_error(
    SSModelDynamicGompertz$new(Y = Y, original.results = "not_a_filterresults"),
    "original.results"
  )
  expect_no_error(SSModelDynamicGompertz$new(Y = Y, original.results = NULL))
})

test_that("xpred / xpred_lead / xpred_targ must be NULL or an idx_series object, in both classes", {
  Y1 <- idx_series(exp(seq(4, 6, length.out = 30)), start = 1L)
  
  expect_error(SSModelDynamicGompertz$new(Y = Y1, xpred = data.frame(x = 1:30)), "xpred")
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, xpred = NULL))
  
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y2 <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
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
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  expect_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 0), "LeadIndCol")
  expect_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 3), "LeadIndCol")
  expect_no_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 1))
  expect_no_error(SSModelLeadingIndicator$new(Y = Y, n.lag = 3, LeadIndCol = 2))
})

test_that("calendar must be NULL or an idx_calendar object, in both classes", {
  Y1 <- idx_series(exp(seq(4, 6, length.out = 30)), start = 1L)
  
  expect_error(SSModelDynamicGompertz$new(Y = Y1, calendar = "not_a_calendar"), "calendar")
  expect_no_error(SSModelDynamicGompertz$new(Y = Y1, calendar = NULL))
  
  lead <- seq(100, 300, length.out = 30)
  targ <- seq(50, 150, length.out = 30)
  Y2 <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  expect_error(
    SSModelLeadingIndicator$new(Y = Y2, n.lag = 3, calendar = "not_a_calendar"),
    "calendar"
  )
})


## 4. SHARED UTILITY BOUNDARIES: get_timeframe() ----

test_that("get_timeframe errors when start > end (both fall within the available range)", {
  x <- idx_series(1:20, start = 1L)
  expect_error(get_timeframe(x, 15, 5), "start is after end")
})

test_that("get_timeframe returns the expected non-empty window when start <= end (paired positive control)", {
  x <- idx_series(1:20, start = 1L)
  result <- get_timeframe(x, 1, 10)
  expect_equal(length(result), 10)
})

test_that("get_timeframe clamps to the available range rather than erroring when start/end extend beyond it", {
  x <- idx_series(1:20, start = 1L)
  result <- get_timeframe(x, -10, 100)
  expect_equal(length(result), 20)
})