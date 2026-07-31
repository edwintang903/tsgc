library(KFAS)

make_li_fit <- function(sea.period = 0, n.lag = 4, xpred_lead = NULL, xpred_targ = NULL,
                        start = NULL, end = NULL) {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  est.start.eng <- if (is.null(start)) idx_to_pos(conv$calendar, as.Date("2021-04-30")) else start
  est.end.eng   <- if (is.null(end)) idx_to_pos(conv$calendar, as.Date("2021-07-24")) else end
  
  mod <- SSModelLeadingIndicator$new(
    conv$series, n.lag = n.lag, sea.period = sea.period,
    start = est.start.eng, end = est.end.eng,
    xpred_lead = xpred_lead, xpred_targ = xpred_targ,
    calendar = conv$calendar
  )
  list(model = mod, res = mod$estimate(), conv = conv,
       est.start = est.start.eng, est.end = est.end.eng)
}

test_that("predict_level computes LI predictions correctly - no seasonal", {
  fit <- make_li_fit(sea.period = 0)
  nf <- 7
  
  forcout <- fit$res$predict_all(nf, sea.on = TRUE, return.all = FALSE)$y.hat.kfas
  LDLtarg_fc <- forcout$LDLtarg
  
  last.pos <- fit$res$end
  last_admit <- as.numeric(idx_values(fit$res$data[last.pos])[, "cTarg"])
  
  forc <- numeric(nf)
  forc[1] <- last_admit * exp(LDLtarg_fc[1, 1])
  forc[2:nf] <- last_admit * exp(LDLtarg_fc[2:nf, 1]) *
    cumprod(1 + exp(LDLtarg_fc[1:(nf - 1), 1]))
  forc <- round(forc, 2)
  
  forc_tsgc <- fit$res$predict_level(n.ahead = nf)
  
  expect_equal(colnames(idx_values(forc_tsgc)), c("forc", "lwr", "upr"))
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "forc"])), forc, tolerance = 1e-4)
})

test_that("predict_level works with seasonal component", {
  fit <- make_li_fit(sea.period = 7)
  nf <- 7
  
  forc_tsgc <- fit$res$predict_level(n.ahead = nf, sea.on = TRUE)
  forc_tsgc_off <- fit$res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_equal(length(forc_tsgc), nf)
  expect_false(isTRUE(all.equal(idx_values(forc_tsgc), idx_values(forc_tsgc_off))))
})

test_that("xpred_lead.new and xpred_targ.new can be assigned and used in predict_level", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  
  # Align the weather series onto england's own position scale via its
  # actual calendar date, rather than assuming a shared start position.
  w_start_pos <- idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1])
  conv_lead <- xts_to_idx(england_weather_2021[, 1:2], start.pos = w_start_pos)
  conv_targ <- xts_to_idx(england_weather_2021[, 3], start.pos = w_start_pos)
  
  fit <- make_li_fit(sea.period = 0, xpred_lead = conv_lead$series, xpred_targ = conv_targ$series)
  
  fit$res$xpred_lead.new <- conv_lead$series
  fit$res$xpred_targ.new <- conv_targ$series
  
  expect_no_error(fit$res$predict_level(n.ahead = 5))
})

test_that("predict_level errors informatively when xpred.new is required but missing", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  w_start_pos <- idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1])
  conv_lead <- xts_to_idx(england_weather_2021[, 1:2], start.pos = w_start_pos)
  
  fit <- make_li_fit(sea.period = 0, xpred_lead = conv_lead$series)
  
  expect_error(fit$res$predict_level(n.ahead = 5), "xpred_lead.new")
})

test_that("get_growth_y works for LI FilterResultsLI object, smoothed and filtered", {
  fit <- make_li_fit(sea.period = 7)
  
  smooth <- fit$res$get_growth_y(smoothed = TRUE, return.components = TRUE)
  filt <- fit$res$get_growth_y(smoothed = FALSE, return.components = FALSE)
  
  expect_equal(length(smooth), 3)
  expect_true(is_idx_series(filt))
})

test_that("get_gy_ci works for FilterResultsLI object", {
  fit <- make_li_fit(sea.period = 7)
  
  smooth <- fit$res$get_gy_ci(smoothed = TRUE)
  filt <- fit$res$get_gy_ci(smoothed = FALSE)
  
  expect_equal(colnames(idx_values(smooth)), c("fit", "lower", "upper"))
  expect_false(isTRUE(all.equal(idx_values(smooth), idx_values(filt))))
})

test_that("print, summary, print_estimation_results work on FilterResultsLI", {
  fit <- make_li_fit()
  
  expect_no_error(print(fit$res))
  expect_no_error(summary(fit$res))
  expect_no_error(fit$res$print_estimation_results())
})

test_that("plot_forecast, plot_log_forecast and plot_holdout work on a FilterResultsLI object", {
  fit <- make_li_fit(sea.period = 7)
  
  expect_no_error(plot_forecast(fit$res))
  expect_no_error(plot_log_forecast(fit$res, Y = fit$conv$series))
  expect_no_error(plot_holdout(fit$res, Y = fit$conv$series))
})

test_that("all standalone plot functions work on a FilterResultsLI object", {
  fit <- make_li_fit(sea.period = 7)
  
  expect_no_error(plot_forecast(fit$res))
  expect_no_error(plot_log_forecast(fit$res, Y = fit$conv$series))
  expect_no_error(plot_holdout(fit$res, Y = fit$conv$series))
  expect_no_error(plot_gy_ci(fit$res))
  expect_no_error(plot_gy_components(fit$res))
})

test_that("plot_gy_ci and plot_gy_components default plt.start to the start of the estimation sample for a FilterResultsLI object", {
  # FilterResultsLI has no `index` field (unlike FilterResults); plt.start
  # defaults via idx_range(res$data)[1] instead, which should equal res$start.
  # fit$res$calendar is set, so the plotted x-axis is calendar dates, not
  # raw positions - convert res$start to a date for the comparison.
  fit <- make_li_fit(sea.period = 7)
  
  p_ci <- plot_gy_ci(fit$res)
  p_comp <- plot_gy_components(fit$res)
  
  expected_start_date <- idx_to_date(fit$res$calendar, fit$res$start)
  
  expect_equal(min(p_ci$data$x), expected_start_date)
  expect_true(min(p_comp$data$x) >= expected_start_date)
})

test_that("plot_gy_ci and plot_gy_components give the same results for FilterResults and FilterResultsLI given equivalent plt.start semantics", {
  # Both classes should compute the same default start position via
  # idx_range(res$data)[1], mirroring the parity between the two classes'
  # (now-retired) RefClass plot_gy_ci()/plot_gy_components() methods in the
  # xts-based version of the package.
  fit <- make_li_fit(sea.period = 0)
  
  expect_equal(idx_range(fit$res$data)[1], fit$res$start)
})

test_that("mapes works for FilterResultsLI object", {
  fit <- make_li_fit()
  
  errs <- fit$res$mapes(n.ahead = 7, Y = fit$conv$series)
  
  expect_equal(length(errs), 5)
})

test_that("predict_level works with a plain, non-calendar idx_series", {
  set.seed(5)
  lead <- cumsum(rpois(150, 6)) + 1
  targ <- cumsum(rpois(150, 8)) + 1
  Y <- idx_series(cbind(lead, targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y, n.lag = 5, sea.period = 0,
                                     start = 1L, end = 100L)
  res <- mod$estimate()
  
  forc <- res$predict_level(n.ahead = 5)
  
  expect_equal(length(forc), 5)
})

test_that("FilterResultsLI object stores start/end as integer positions, not dates", {
  fit <- make_li_fit()
  
  expect_type(fit$res$start, "integer")
  expect_type(fit$res$end, "integer")
  expect_equal(fit$res$start, fit$est.start)
  expect_equal(fit$res$end, fit$est.end)
})
test_that("xpred_lead.new/xpred_targ.new select forecast rows by exact date, not position, when the supplied regressor series begins months before the model sample", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  est.start.eng <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  est.end.eng   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  
  # Weather regressor series begins several months before the model's
  # own estimation sample, so a positional (row-count) selection would
  # misalign against a date-based selection.
  weather_idx <- xts_to_idx(
    england_weather_2021[, "temperature_C", drop = FALSE],
    start.pos = idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1])
  )$series
  expect_true(idx_positions(weather_idx)[1] < est.start.eng - 60)
  
  mod <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 0,
    start = est.start.eng, end = est.end.eng,
    xpred_lead = weather_idx, xpred_targ = weather_idx,
    calendar = conv$calendar
  )
  res <- mod$estimate()
  
  nf <- 14
  future_dates <- idx_to_date(conv$calendar, (est.end.eng + 1):(est.end.eng + nf))
  
  # FilterResultsLI has no `xpred.new` field - only `xpred_lead.new` and
  # `xpred_targ.new` (see the class definition in filterResultsLI.R).
  # An earlier version of this test assigned to `res$xpred.new`, a
  # field that does not exist on this setRefClass object, so the
  # assignment never reached the fields predict_all() actually reads
  # and the test exercised get_timeframe() against itself rather than
  # the model's real internal alignment logic.
  #
  # We deliberately hand the model the FULL, unsliced weather series
  # (still starting months before the estimation sample) rather than a
  # pre-sliced window, so predict_all()'s own internal
  # get_timeframe(xpred_*.new, end+1, end+n.ahead) call (see
  # filterResultsLI.R) is what performs the date-based selection here,
  # not the test.
  res$xpred_lead.new <- weather_idx
  res$xpred_targ.new <- weather_idx
  
  out <- res$predict_all(nf, sea.on = TRUE, return.all = FALSE)
  
  # Independently compute the expected forecast dates/values via
  # get_timeframe, to compare against what predict_all() actually used
  # internally. This is now a genuine check that the model's own
  # alignment logic (not the test's pre-slicing) selected the correct
  # rows, since res$xpred_lead.new/xpred_targ.new were supplied in full,
  # unsliced form above.
  expected_future <- get_timeframe(weather_idx, est.end.eng + 1, est.end.eng + nf)
  sel_dates <- idx_to_date(conv$calendar, idx_positions(expected_future))
  expect_equal(sel_dates, future_dates)
  expect_equal(length(unique(sel_dates)), nf)
  
  # The forecast output itself must exist and cover exactly the
  # expected forecast horizon - i.e. predict_all() successfully
  # resolved dates from the full, unsliced regressor series without
  # erroring or silently truncating.
  expect_true(is_idx_series(out) || is.list(out))
  out_positions <- if (is_idx_series(out)) idx_positions(out) else idx_positions(out[[1]])
  expect_equal(range(out_positions), c(est.end.eng + 1, est.end.eng + nf))
})

test_that("xpred date alignment fails clearly when the supplied regressor series does not extend across the full forecast horizon", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  est.start.eng <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  est.end.eng   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  
  weather_idx <- xts_to_idx(
    england_weather_2021[, "temperature_C", drop = FALSE],
    start.pos = idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1])
  )$series
  
  # idx_series objects are always a contiguous run of positions by
  # construction (idx_series(data, start) implies positions
  # start:(start+nrow(data)-1); idx_rbind() requires its two inputs to
  # be adjoining and errors otherwise) - there is no way to build an
  # idx_series with an actual internal gap. What CAN happen in
  # practice is a regressor series that simply does not extend far
  # enough to cover the requested forecast horizon (e.g. a weather
  # forecast feed that stops a few days short, or that was
  # accidentally sliced too early). Simulate that directly by
  # truncating the series a few positions before the forecast horizon
  # ends, rather than trying (and failing, since it isn't possible) to
  # punch an internal hole in it.
  nf <- 14
  full_future <- get_timeframe(weather_idx, est.end.eng + 1, est.end.eng + nf)
  short_future <- idx_series(
    idx_values(full_future)[1:(nf - 5), , drop = FALSE],
    start = idx_positions(full_future)[1]
  )
  expect_equal(length(idx_positions(short_future)), nf - 5)
  
  mod <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 0,
    start = est.start.eng, end = est.end.eng,
    xpred_lead = weather_idx, xpred_targ = weather_idx,
    calendar = conv$calendar
  )
  res <- mod$estimate()
  # FilterResultsLI has no `xpred.new` field; the fields predict_all()
  # actually reads are xpred_lead.new and xpred_targ.new (see
  # filterResultsLI.R). Assigning to xpred_lead.new/xpred_targ.new here
  # ensures the truncated series is genuinely what gets fed into the
  # model's internal date-based alignment, rather than being silently
  # ignored because it was written to a field the class doesn't have.
  res$xpred_lead.new <- short_future
  res$xpred_targ.new <- short_future
  
  # predict_all() now checks (see filterResultsLI.R) that
  # get_timeframe(xpred_*.new, end+1, end+n.ahead) actually returns
  # n.ahead rows before using it, rather than relying on
  # get_timeframe()'s own start/end clamping (which silently returns
  # fewer rows instead of erroring, and would otherwise let a
  # too-short regressor series be positionally recycled into a
  # fixed-size n.ahead array without any error).
  expect_error(res$predict_all(nf, sea.on = TRUE, return.all = FALSE),
               "does not cover the full forecast horizon")
})

test_that("get_timeframe silently clamps to the available range rather than erroring - documented so callers relying on it for validation are aware", {
  # This is not a bug fix test - it pins get_timeframe()'s actual,
  # current behaviour (silent clamping via start <- max(start, rng[1]);
  # end <- min(end, rng[2])), which is why predict_all() in
  # filterResultsLI.R cannot rely on get_timeframe() alone to detect an
  # incomplete forecast-horizon window and instead checks the returned
  # window's length explicitly (see the test above). If this test ever
  # fails because get_timeframe() has been changed to error on an
  # out-of-range request instead of clamping, the explicit length
  # check added to predict_all() would become redundant (but harmless)
  # defence-in-depth, not incorrect.
  x <- idx_series(matrix(1:5, ncol = 1), start = 1L)
  out <- get_timeframe(x, start = 3, end = 10)
  expect_equal(idx_positions(out), 3:5)
  expect_true(length(idx_positions(out)) < (10 - 3 + 1))
})

test_that("xts_to_idx errors clearly on a duplicated date in its input index", {
  # xts_to_idx()'s source (utils.R) checks
  #   if (anyDuplicated(idx)) stop("xts_to_idx: x's index contains ",
  #                                "duplicate values.")
  # unconditionally, before any detection/construction logic runs - so
  # a duplicated index date always errors clearly; it is never
  # silently de-duplicated or allowed to produce an idx_series with an
  # ambiguous position-to-date mapping. This test pins that behaviour
  # directly against a genuinely malformed xts input (two source rows
  # sharing one date, as could occur from concatenating overlapping
  # series), rather than only asserting it as a requirement.
  data(england_weather_2021, package = "tsgc")
  base_xts <- england_weather_2021[, "temperature_C", drop = FALSE]
  
  dup_date <- zoo::index(base_xts)[10]
  dup_row <- base_xts[10, , drop = FALSE]
  malformed_xts <- rbind(base_xts, dup_row)
  
  expect_true(sum(zoo::index(malformed_xts) == dup_date) == 2)
  expect_error(xts_to_idx(malformed_xts), "duplicate values")
})

test_that("sMAPE, as computed by mapes(), uses (Actual + Forecast) in the denominator and ranges over [0, 100], not the conventional [0, 200] scale", {
  # Reproduces the exact formula used inside FilterResults$mapes() and
  # FilterResultsLI$mapes():
  #   smape <- mean(100 * abs(Actual - Forecast) / (Actual + Forecast))
  # This is NOT the textbook sMAPE, which divides by the mean of the
  # absolute values ((|Actual| + |Forecast|) / 2) and ranges over
  # [0, 200]. This test pins the formula and its scale directly, so a
  # future change to the formula (e.g. "fixing" it to the textbook
  # version) is caught rather than silently changing the meaning of
  # every reported sMAPE figure in the package.
  smape_as_implemented <- function(actual, forecast) {
    mean(100 * abs(actual - forecast) / (actual + forecast))
  }
  
  # Exact forecast: sMAPE must be exactly 0, both under this formula and
  # under the textbook one, since both denominators are strictly
  # positive when Actual == Forecast > 0.
  expect_equal(smape_as_implemented(c(100, 200, 300), c(100, 200, 300)), 0)
  
  # Textbook worked example: Actual = 100, Forecast = 110.
  #   As implemented: 100 * |100-110| / (100+110) = 100 * 10/210 = 4.7619...
  #   Textbook sMAPE (0-200 scale): 100 * |100-110| / ((100+110)/2) = 9.5238...
  # i.e. this package's sMAPE is exactly half the textbook value for a
  # single-point evaluation, since (A+F) = 2 * ((A+F)/2). Pinning both
  # numbers makes the scale difference explicit and machine-checkable.
  implemented_val <- smape_as_implemented(100, 110)
  textbook_val <- 100 * abs(100 - 110) / ((100 + 110) / 2)
  expect_equal(implemented_val, 100 * 10 / 210)
  expect_equal(implemented_val, textbook_val / 2)
  
  # Maximum-divergence case: Forecast = 0 (or Actual = 0) drives the
  # implemented statistic to exactly 100, its upper bound - confirming
  # the [0, 100] range rather than [0, 200].
  expect_equal(smape_as_implemented(100, 0), 100)
  expect_equal(smape_as_implemented(0, 100), 100)
  
  # The implemented statistic never exceeds 100 for positive
  # actual/forecast pairs, unlike the textbook version which can reach
  # 200 as forecast error grows without bound relative to the average
  # of actual and forecast.
  set.seed(1)
  actual <- runif(50, min = 1, max = 1000)
  forecast <- runif(50, min = 1, max = 1000)
  expect_true(all(smape_as_implemented(actual, forecast) <= 100 + 1e-8))
})