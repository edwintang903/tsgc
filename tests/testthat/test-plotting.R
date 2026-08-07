test_that("forecast plots keep angled date labels clear of the axis title", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  model <- tsgc::SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49), q = 0.005, calendar = conv$calendar)
  res <- model$estimate()
  
  plots <- list(
    plot_forecast(res),
    plot_log_forecast(res, conv$series)
  )
  for (p in plots) {
    expect_equal(p$theme$axis.text.x$angle, 45)
    expect_equal(p$theme$axis.text.x$hjust, 1)
    expect_equal(p$theme$axis.text.x$vjust, 1)
    expect_equal(as.numeric(p$theme$axis.text.x$margin[1]), 6)
    expect_equal(as.numeric(p$theme$axis.title.x$margin[1]), 18)
  }
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

test_that("plot_holdout keeps angled date labels clear of the axis title", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = conv$calendar)
  res <- model$estimate()
  # Plot forecasts and outcomes over evaluation period
  p <- plot_holdout(res = res, Y = conv$series)
  expect_equal(p$theme$axis.text.x$angle, 45)
  expect_equal(p$theme$axis.text.x$hjust, 1)
  expect_equal(p$theme$axis.text.x$vjust, 1)
  expect_equal(as.numeric(p$theme$axis.text.x$margin[1]), 6)
  expect_equal(as.numeric(p$theme$axis.title.x$margin[1]), 18)
})

test_that("plot_compare_forecast uses the forecast date-label style", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  
  model1 <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = conv$calendar)
  model2 <- SSModelDynamicGompertz$new(Y = conv$series, q = NULL, end = end.pos, calendar = conv$calendar)
  res1 <- estimate(model1)
  res2 <- estimate(model2)
  
  p <- plot_compare_forecast(
    list(q_fixed = res1, q_estimated = res2), n.ahead = 7
  )
  expect_equal(p$theme$axis.text.x$angle, 45)
  expect_equal(p$theme$axis.text.x$hjust, 1)
  expect_equal(p$theme$axis.text.x$vjust, 1)
  expect_equal(as.numeric(p$theme$axis.text.x$margin[1]), 6)
  expect_equal(as.numeric(p$theme$axis.title.x$margin[1]), 18)
})

test_that("plot.SSModelLeadingIndicator uses the forecast date-label style", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  model <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 7,
    start = idx_to_pos(conv$calendar, as.Date("2021-04-30")),
    end = idx_to_pos(conv$calendar, as.Date("2021-07-24")),
    calendar = conv$calendar
  )
  
  p <- plot(model)
  expect_equal(p$theme$axis.text.x$angle, 45)
  expect_equal(p$theme$axis.text.x$hjust, 1)
  expect_equal(p$theme$axis.text.x$vjust, 1)
  expect_equal(as.numeric(p$theme$axis.text.x$margin[1]), 6)
  expect_equal(as.numeric(p$theme$axis.title.x$margin[1]), 18)
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
test_that("plot_r0() works on a fitted SSModelDynamicGompertz result", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2021-05-03"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
                                      end = end.pos, calendar = conv$calendar)
  res <- estimate(model)
  
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 7))
})

test_that("plot_r0() works on a fitted SSModelLeadingIndicator result", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  start.pos <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  end.pos   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  
  model <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 7,
    start = start.pos, end = end.pos, calendar = conv$calendar
  )
  res <- model$estimate()
  
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 7))
})

test_that("plot_r0() falls back to plain integer x-axis positions when no calendar is supplied", {
  set.seed(1)
  Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1L)
  model <- SSModelDynamicGompertz$new(Y = Y, q = NULL, end = 100)
  res <- estimate(model)
  
  expect_no_error(plot_r0(res, gen_int = 5, n.ahead = 7))
})

test_that("plot_r0() errors when res is not a FilterResults or FilterResultsLI object", {
  expect_error(plot_r0(list(), gen_int = 5), "FilterResults")
})

test_that("plot_r0() respects axis modes and honours n.ahead/smoothed/confidence.level", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2021-05-03"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
                                      end = end.pos, calendar = conv$calendar)
  res <- estimate(model)
  
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 5,
                          axis = idx_axis_opts(mode = "position")))
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 5,
                          axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 5,
                          smoothed = TRUE, confidence.level = 0.9))
})

## ----------------------------------------------------------------------
## Non-posixct calendar coverage: every plot.*/plot_* entry point
## ----------------------------------------------------------------------

test_that("plot(model) methods work with a non-posixct calendar (Gompertz and LeadingIndicator)", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  model <- SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49),
                                      q = 0.005, calendar = cal_np)
  expect_no_error(plot(model))
  expect_no_error(plot(model, axis = idx_axis_opts(mode = "position")))
  expect_no_error(plot(model, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot(model, axis = idx_axis_opts(mode = "time_since")))
  
  data(england, package = "tsgc")
  conv2 <- xts_to_idx(england[, 1:2])
  cal_np2 <- idx_calendar(anchor = zoo::index(england)[1], anchor_pos = 1L,
                          amount = 1, unit = "days", posixct = FALSE)
  model_li <- SSModelLeadingIndicator$new(
    conv2$series, n.lag = 4, sea.period = 7,
    start = idx_to_pos(conv2$calendar, as.Date("2021-04-30")),
    end = idx_to_pos(conv2$calendar, as.Date("2021-07-24")),
    calendar = cal_np2
  )
  expect_no_error(plot(model_li))
  expect_no_error(plot(model_li, axis = idx_axis_opts(mode = "steps")))
})

test_that("plot_forecast() and plot_log_forecast() work with a non-posixct calendar, all axis modes", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE,
                         anchor_name = "first recorded case")
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = cal_np)
  res <- estimate(model)
  
  expect_no_error(plot_forecast(res))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "position")))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "time_since")))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "steps", info_box = TRUE)))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "time_since", info_box = TRUE, pattern_n = 1)))
  # "auto" with a non-posixct calendar must resolve to "position", not "date".
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "auto")))
  
  expect_no_error(plot_log_forecast(res, conv$series))
  expect_no_error(plot_log_forecast(res, conv$series, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_log_forecast(res, conv$series, axis = idx_axis_opts(mode = "time_since")))
})

test_that("plot_forecast() and plot_log_forecast() work with a non-posixct calendar for SSModelLeadingIndicator", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  cal_np <- idx_calendar(anchor = zoo::index(england)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  start.pos <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  end.pos   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  model <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 7,
    start = start.pos, end = end.pos, calendar = cal_np
  )
  res <- estimate(model)
  
  expect_no_error(plot_forecast(res))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "time_since")))
  expect_no_error(plot_log_forecast(res, conv$series))
  expect_no_error(plot_log_forecast(res, conv$series, axis = idx_axis_opts(mode = "time_since")))
})

test_that("plot_gy_components() and plot_gy_ci() work with a non-posixct calendar, all axis modes", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  model <- SSModelDynamicGompertz$new(Y = get_timeframe(conv$series, conv$series$start, conv$series$start + 49),
                                      q = 0.005, calendar = cal_np)
  res <- estimate(model)
  
  expect_no_error(plot_gy_components(res))
  expect_no_error(plot_gy_components(res, axis = idx_axis_opts(mode = "position")))
  expect_no_error(plot_gy_components(res, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_gy_components(res, axis = idx_axis_opts(mode = "time_since")))
  expect_no_error(plot_gy_components(res, smoothed = TRUE))
  
  expect_no_error(plot_gy_ci(res))
  expect_no_error(plot_gy_ci(res, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_gy_ci(res, axis = idx_axis_opts(mode = "time_since")))
  expect_no_error(plot_gy_ci(res, smoothed = TRUE))
  
  data(england, package = "tsgc")
  conv2 <- xts_to_idx(england[, 1:2])
  cal_np2 <- idx_calendar(anchor = zoo::index(england)[1], anchor_pos = 1L,
                          amount = 1, unit = "days", posixct = FALSE)
  model_li <- SSModelLeadingIndicator$new(
    conv2$series, n.lag = 4, sea.period = 7,
    start = idx_to_pos(conv2$calendar, as.Date("2021-04-30")),
    end = idx_to_pos(conv2$calendar, as.Date("2021-07-24")),
    calendar = cal_np2
  )
  res_li <- estimate(model_li)
  expect_no_error(plot_gy_ci(res_li))
  expect_no_error(plot_gy_components(res_li, axis = idx_axis_opts(mode = "time_since")))
})

test_that("plot_holdout() works with a non-posixct calendar for both model types", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = cal_np)
  res <- estimate(model)
  
  expect_no_error(plot_holdout(res = res, Y = conv$series))
  expect_no_error(plot_holdout(res = res, Y = conv$series, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_holdout(res = res, Y = conv$series, axis = idx_axis_opts(mode = "time_since")))
  
  data(england, package = "tsgc")
  conv2 <- xts_to_idx(england[, 1:2])
  cal_np2 <- idx_calendar(anchor = zoo::index(england)[1], anchor_pos = 1L,
                          amount = 1, unit = "days", posixct = FALSE)
  start.pos2 <- idx_to_pos(conv2$calendar, as.Date("2021-04-30"))
  end.pos2   <- idx_to_pos(conv2$calendar, as.Date("2021-07-24"))
  model_li <- SSModelLeadingIndicator$new(
    conv2$series, n.lag = 4, sea.period = 7,
    start = start.pos2, end = end.pos2, calendar = cal_np2
  )
  res_li <- estimate(model_li)
  expect_no_error(plot_holdout(res = res_li, Y = conv2$series))
  expect_no_error(plot_holdout(res = res_li, Y = conv2$series, axis = idx_axis_opts(mode = "steps")))
})

test_that("plot_compare_forecast() works across two models sharing a non-posixct calendar", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  
  model1 <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = cal_np)
  model2 <- SSModelDynamicGompertz$new(Y = conv$series, q = NULL, end = end.pos, calendar = cal_np)
  res1 <- estimate(model1)
  res2 <- estimate(model2)
  
  expect_no_error(plot_compare_forecast(list(q_fixed = res1, q_estimated = res2), n.ahead = 7))
  expect_no_error(plot_compare_forecast(list(q_fixed = res1, q_estimated = res2), n.ahead = 7,
                                        axis = idx_axis_opts(mode = "steps")))
})

test_that("plot_r0() works with a non-posixct calendar for both model types", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2021-05-03"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = cal_np)
  res <- estimate(model)
  
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 7))
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 7, axis = idx_axis_opts(mode = "position")))
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 7, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 7, axis = idx_axis_opts(mode = "time_since")))
  expect_no_error(plot_r0(res, gen_int = 4, n.ahead = 5, smoothed = TRUE, confidence.level = 0.9))
  
  data(england, package = "tsgc")
  conv2 <- xts_to_idx(england[, 1:2])
  cal_np2 <- idx_calendar(anchor = zoo::index(england)[1], anchor_pos = 1L,
                          amount = 1, unit = "days", posixct = FALSE)
  model_li <- SSModelLeadingIndicator$new(
    conv2$series, n.lag = 4, sea.period = 7,
    start = idx_to_pos(conv2$calendar, as.Date("2021-04-30")),
    end = idx_to_pos(conv2$calendar, as.Date("2021-07-24")),
    calendar = cal_np2
  )
  res_li <- estimate(model_li)
  expect_no_error(plot_r0(res_li, gen_int = 4, n.ahead = 7))
  expect_no_error(plot_r0(res_li, gen_int = 4, n.ahead = 7, axis = idx_axis_opts(mode = "steps")))
})

test_that("plot_* functions work with a fully non-calendar (numeric-anchor) idx_calendar", {
  set.seed(7)
  Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1L)
  cal_ps <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                         unit = "picoseconds", posixct = FALSE)
  model <- SSModelDynamicGompertz$new(Y = Y, q = NULL, end = 100, calendar = cal_ps)
  res <- estimate(model)
  
  expect_no_error(plot_forecast(res))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "steps")))
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "time_since")))
  expect_no_error(plot_gy_ci(res))
  expect_no_error(plot_gy_components(res))
  expect_no_error(plot_r0(res, gen_int = 5, n.ahead = 7))
  expect_no_error(plot_holdout(res = res, Y = Y))
})

test_that("plot_forecast's x-axis is plain numeric positions by default under posixct = FALSE (mode='auto')", {
  set.seed(24)
  Y <- idx_series(cumsum(rpois(60, 8)) + 1, start = 1L)
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 50, calendar = cal_np)
  res <- estimate(model)
  
  p <- plot_forecast(res)
  expect_equal(p$labels$x, "Position")
  expect_true(is.numeric(p$data$x))
})