library(KFAS)

test_that("forecasting does not truncate supplied future regressors", {
  set.seed(3)
  Y <- idx_series(cumsum(rpois(130, 8)) + 1)
  xpred <- idx_series(matrix(rnorm(130), ncol = 1))
  model <- SSModelDynamicGompertz(
    Y, q = 0.005, sea.period = 0, start = 1, end = 100, xpred = xpred
  )
  res <- estimate(model)
  res$xpred.new <- xpred
  original_xpred <- res$xpred.new

  expect_no_error(res$predict_all(5))
  expect_identical(res$xpred.new, original_xpred)
  expect_no_error(res$predict_all(10))
  expect_identical(res$xpred.new, original_xpred)
})

test_that("predict_level computes predictions correctly - no seasonal", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 0,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf,
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  cp <- cumprod(1 + exp(delta_fit))
  mult <- c(1, cp[1:(nf - 1)])
  forc <- rep(YT, nf) * exp(delta_fit) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  cp_lwr <- cumprod(1 + exp(delta_lwr))
  mult_lwr <- c(1, cp_lwr[1:(nf - 1)])
  forc_lwr <- rep(YT, nf) * exp(delta_lwr) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions of cumulated variable correctly - no seasonal", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 0,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf,
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  mult <- cumprod(1 + exp(delta_fit))
  forc <- rep(YT, nf) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  mult_lwr <- cumprod(1 + exp(delta_lwr))
  forc_lwr <- rep(YT, nf) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, return.diff = FALSE)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 1)
})

test_that("predict_level computes predictions correctly - seasonal, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf,
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  cp <- cumprod(1 + exp(delta_fit))
  mult <- c(1, cp[1:(nf - 1)])
  forc <- rep(YT, nf) * exp(delta_fit) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  cp_lwr <- cumprod(1 + exp(delta_lwr))
  mult_lwr <- c(1, cp_lwr[1:(nf - 1)])
  forc_lwr <- rep(YT, nf) * exp(delta_lwr) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = TRUE)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal but sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf,
                        interval = c("confidence"), level = 0.68,
                        states = c("level"))
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  cp <- cumprod(1 + exp(delta_fit))
  mult <- c(1, cp[1:(nf - 1)])
  forc <- rep(YT, nf) * exp(delta_fit) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  cp_lwr <- cumprod(1 + exp(delta_lwr))
  mult_lwr <- c(1, cp_lwr[1:(nf - 1)])
  forc_lwr <- rep(YT, nf) * exp(delta_lwr) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + AR1, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf,
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  cp <- cumprod(1 + exp(delta_fit))
  mult <- c(1, cp[1:(nf - 1)])
  forc <- rep(YT, nf) * exp(delta_fit) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  cp_lwr <- cumprod(1 + exp(delta_lwr))
  mult_lwr <- c(1, cp_lwr[1:(nf - 1)])
  forc_lwr <- rep(YT, nf) * exp(delta_lwr) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = TRUE)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + AR1, sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf,
                        interval = c("confidence"), level = 0.68,
                        states = c("level", "custom"))
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  cp <- cumprod(1 + exp(delta_fit))
  mult <- c(1, cp[1:(nf - 1)])
  forc <- rep(YT, nf) * exp(delta_fit) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  cp_lwr <- cumprod(1 + exp(delta_lwr))
  mult_lwr <- c(1, cp_lwr[1:(nf - 1)])
  forc_lwr <- rep(YT, nf) * exp(delta_lwr) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + xpred + AR1, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  w_start_pos <- idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = w_start_pos)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, xpred = conv_w$series, sea.period = 7,
    start = est.start, end = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  res$xpred.new <- conv_w$series
  
  f.start <- est.end + 1
  f.end <- est.end + nf
  
  new_weather <- idx_values(get_timeframe(conv_w$series, f.start, f.end))
  
  Qt.slope <- res$output$model$Q[2, 2, 1]
  Qt.seas <- res$output$model$Q[3, 3, 1]
  Qt.ar1 <- res$output$model$Q[9, 9, 1]
  Ht <- res$output$model$H[1, 1, 1]
  
  new_model <- SSModel(formula = matrix(rep(NA, nf), ncol = 1) ~
                         SSMtrend(degree = 2, Q = list(matrix(0),
                                                       matrix(Qt.slope)))
                       + SSMseasonal(period = 7, Q = Qt.seas,
                                     sea.type = "trigonometric")
                       + SSMregression(~new_weather)
                       + SSMcustom(Z = 1, T = 1, R = 1, Q = Qt.ar1,
                                   state_names = "ar1"),
                       H = matrix(Ht))
  
  delta_pred <- predict(res$output$model, newdata = new_model,
                        interval = c("confidence"), level = 0.68, states = c("all"))
  
  delta_fit <- as.vector(delta_pred[, "fit"])
  last.pos <- idx_range(model$Y)[2]
  YT <- as.numeric(idx_values(model$Y[last.pos]))
  cp <- cumprod(1 + exp(delta_fit))
  mult <- c(1, cp[1:(nf - 1)])
  forc <- rep(YT, nf) * exp(delta_fit) * mult
  
  delta_lwr <- as.vector(delta_pred[, "lwr"])
  cp_lwr <- cumprod(1 + exp(delta_lwr))
  mult_lwr <- c(1, cp_lwr[1:(nf - 1)])
  forc_lwr <- rep(YT, nf) * exp(delta_lwr) * mult_lwr
  
  delta_upr <- as.vector(delta_pred[, "upr"])
  cp_upr <- cumprod(1 + exp(delta_upr))
  mult_upr <- c(1, cp_upr[1:(nf - 1)])
  forc_upr <- rep(YT, nf) * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = TRUE)
  
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "fit"])), forc)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "lower"])), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(idx_values(forc_tsgc)[, "upper"])), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + xpred + AR1, sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  w_start_pos <- idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = w_start_pos)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, xpred = conv_w$series, sea.period = 7,
    start = est.start, end = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  res$xpred.new <- conv_w$series
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_true(is_idx_series(forc_tsgc))
  expect_equal(length(forc_tsgc), nf)
  expect_equal(colnames(idx_values(forc_tsgc)), c("fit", "lower", "upper"))
})

test_that("predict_level works - quarterly data", {
  data(nintendo_sales, package = "tsgc")
  # nintendo_sales is indexed by zoo::yearqtr; xts_to_idx() assumes a daily
  # Date/POSIXct index, so build the idx_series/idx_calendar directly with
  # unit = "quarters" instead.
  first_qtr <- zoo::index(nintendo_sales)[1]
  conv <- list(
    series = idx_series(zoo::coredata(nintendo_sales[, 1]), start = 1L),
    calendar = idx_calendar(anchor = zoo::as.Date(first_qtr), anchor_pos = 1L,
                            amount = 1, unit = "quarters", posixct = TRUE)
  )
  
  est.start.q <- idx_to_pos(conv$calendar, as.Date("2006-10-01"))
  est.end.q   <- idx_to_pos(conv$calendar, as.Date("2010-07-01"))
  
  nf <- 4
  
  mod_wii <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, sea.period = 4, start = est.start.q, end = est.end.q
  )
  res_wii <- mod_wii$estimate()
  
  forecasts <- res_wii$predict_level(n.ahead = nf)
  
  expect_equal(length(forecasts), nf)
})

test_that("predict_level works - non-calendar idx_series", {
  set.seed(9)
  Y <- idx_series(cumsum(rpois(60, 8)) + 1, start = 1L)
  
  mod <- tsgc::SSModelDynamicGompertz$new(
    Y = Y, sea.period = 0, start = 1L, end = 50L
  )
  res <- mod$estimate()
  
  forecasts <- res$predict_level(n.ahead = 1)
  
  expect_equal(length(forecasts), 1)
})

test_that("get_growth_y works and toggles smoothed = TRUE, FALSE correctly", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  smooth <- res$get_growth_y(smoothed = TRUE, return.components = TRUE)
  filt <- res$get_growth_y(smoothed = FALSE, return.components = FALSE)
  filt_all <- res$get_growth_y(smoothed = FALSE, return.components = TRUE)
  
  expect_equal(length(smooth), 3)
  expect_true(is_idx_series(filt))
  expect_false(is_idx_series(smooth))
  expect_true(is_idx_series(smooth[[1]]))
  expect_false(isTRUE(all.equal(idx_values(smooth[[1]]), idx_values(filt))))
})

test_that("get_gy_ci works and toggles smoothed = TRUE, FALSE correctly", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  smooth <- res$get_gy_ci(smoothed = TRUE)
  filt <- res$get_gy_ci(smoothed = FALSE)
  
  expect_equal(colnames(idx_values(smooth)), c("fit", "lower", "upper"))
  expect_false(isTRUE(all.equal(idx_values(smooth), idx_values(filt))))
})

test_that("print_estimation_results works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  expect_no_error(res$print_estimation_results())
})

test_that("print and summary methods work", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  expect_no_error(print(res))
  expect_no_error(summary(res))
})

test_that("plot functions work on a FilterResults object (via standalone plotting functions)", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end, calendar = conv$calendar
  )
  res <- estimate(model)
  
  expect_no_error(plot_forecast(res))
  expect_no_error(plot_log_forecast(res, Y = conv$series))
  expect_no_error(plot_gy_ci(res))
  expect_no_error(plot_gy_components(res))
  expect_no_error(plot_holdout(res, Y = conv$series))
})

test_that("xpred.new fails clearly when the supplied regressor series does not extend across the full forecast horizon", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  w_start_pos <- idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = w_start_pos)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  nf <- 14
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, xpred = conv_w$series, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  # get_timeframe() clamps rather than errors on an out-of-range
  # request, so simulate a regressor feed that stops a few days short
  # of the forecast horizon by truncating the available future window,
  # rather than trying to punch an internal gap in an idx_series
  # (which is not possible - idx_series positions are always a
  # contiguous run by construction).
  full_future <- get_timeframe(conv_w$series, est.end + 1, est.end + nf)
  short_future <- idx_series(
    idx_values(full_future)[1:(nf - 5), , drop = FALSE],
    start = idx_positions(full_future)[1]
  )
  expect_equal(length(idx_positions(short_future)), nf - 5)
  
  res$xpred.new <- short_future
  
  # predict_all()/predict_level() now check (see filterResults.R) that
  # get_timeframe(xpred.new, firstpred, firstpred + n.ahead - 1L)
  # actually returns n.ahead rows before using it, rather than relying
  # on get_timeframe()'s own silent start/end clamping, which would
  # otherwise let a too-short regressor series be positionally
  # recycled into the fixed-size n.ahead array with no error.
  expect_error(res$predict_level(n.ahead = nf, sea.on = TRUE),
               "does not cover the full forecast horizon")
})

test_that("mapes works", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005, sea.period = 7,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  errs <- res$mapes(n.ahead = 7, Y = conv$series)
  
  expect_equal(length(errs), 5)
})
