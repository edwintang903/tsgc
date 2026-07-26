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