test_that("Test plot_forecast() works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  model <- tsgc::SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49), q = 0.005, calendar = conv$calendar)
  res <- model$estimate()
  tsgc::plot_forecast(res)
  expect_equal(1, 1)
})

test_that("Test plot_log_forecast() works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  model <- tsgc::SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49), q = 0.005, calendar = conv$calendar)
  res <- model$estimate()
  tsgc::plot_log_forecast(res, conv$series)
  expect_equal(1, 1)
})

test_that("Test plot_gy_components() works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  model <- tsgc::SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49), q = 0.005, calendar = conv$calendar)
  res <- model$estimate()
  tsgc::plot_gy_components(res)
  expect_equal(1, 1)
})

test_that("Test plot_gy_ci() works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  model <- tsgc::SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49), q = 0.005, calendar = conv$calendar)
  res <- model$estimate()
  tsgc::plot_gy_ci(res)
  expect_equal(1, 1)
})

test_that("Test plot_holdout() works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = conv$calendar)
  res <- model$estimate()
  # Plot forecasts and outcomes over evaluation period
  plot_holdout(res = res, Y = conv$series)
  expect_equal(1, 1)
})

test_that("Test plot_compare_forecast() works across two models", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  
  model1 <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = conv$calendar)
  model2 <- SSModelDynamicGompertz$new(Y = conv$series, q = NULL, end = end.pos, calendar = conv$calendar)
  res1 <- estimate(model1)
  res2 <- estimate(model2)
  
  expect_no_error(
    plot_compare_forecast(list(q_fixed = res1, q_estimated = res2), n.ahead = 7)
  )
})

test_that("Plots fall back to plain integer x-axis positions when no calendar is supplied", {
  set.seed(1)
  Y <- idx_series(cumsum(rpois(80, 8)) + 1, start = 1L)
  model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 60)
  res <- estimate(model)
  
  expect_no_error(plot_forecast(res))
  expect_no_error(plot_gy_ci(res))
})

test_that("plot_gy_ci() and plot_gy_components() work for SSModelLeadingIndicator results too", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  start.pos <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  end.pos   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  
  model <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 7,
    start = start.pos, end = end.pos, calendar = conv$calendar
  )
  res <- model$estimate()
  
  expect_no_error(plot_gy_ci(res))
  expect_no_error(plot_gy_components(res))
})