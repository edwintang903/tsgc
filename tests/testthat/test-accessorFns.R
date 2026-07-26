library(tsgc)
library(KFAS)

make_gompertz_fit <- function(end = 100) {
  set.seed(1)
  Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1L)
  model <- SSModelDynamicGompertz$new(
    Y = Y,
    q = 0.005,
    end = end
  )
  estimate(model)
}

test_that("xpred.new field can be assigned directly on a FilterResults object", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  
  conv_g <- xts_to_idx(gauteng)
  
  # Align the weather series onto the SAME position scale as gauteng by
  # converting it to a position using gauteng's own calendar (rather than
  # assuming the two xts objects start on the same calendar date).
  w_start_pos <- idx_to_pos(conv_g$calendar, zoo::index(gauteng_weather_2021)[1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = w_start_pos)
  
  start.pos <- idx_to_pos(conv_g$calendar, as.Date("2021-02-01"))
  end.pos   <- idx_to_pos(conv_g$calendar, as.Date("2021-04-19"))
  
  model_weather <- SSModelDynamicGompertz$new(
    Y = conv_g$series, xpred = conv_w$series,
    start = start.pos, end = end.pos
  )
  res_weather <- estimate(model_weather)
  expect_null(res_weather$xpred.new)
  
  # Feed future weather data into the results object directly via field
  # assignment (the RefClass analogue of the old supply_xpred.new()).
  res_weather$xpred.new <- conv_w$series
  expect_equal(res_weather$xpred.new, conv_w$series)
})

test_that("xpred_lead.new / xpred_targ.new fields can be assigned directly on a FilterResultsLI object", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  
  conv_e <- xts_to_idx(england[, 1:2])
  
  # Align the weather series onto england's own position scale via its
  # actual calendar date, rather than assuming a shared start position.
  w_start_pos <- idx_to_pos(conv_e$calendar, zoo::index(england_weather_2021)[1])
  conv_w <- xts_to_idx(england_weather_2021[, 1:2], start.pos = w_start_pos)
  conv_w3 <- xts_to_idx(england_weather_2021[, 3], start.pos = w_start_pos)
  
  start.pos <- idx_to_pos(conv_e$calendar, as.Date("2021-04-30"))
  end.pos   <- idx_to_pos(conv_e$calendar, as.Date("2021-07-24"))
  
  mod <- SSModelLeadingIndicator$new(
    conv_e$series, n.lag = 4,
    xpred_lead = conv_w$series,
    xpred_targ = conv_w3$series,
    start = start.pos, end = end.pos
  )
  res_lead.x <- estimate(mod)
  
  res_lead.x$xpred_lead.new <- conv_w$series
  res_lead.x$xpred_targ.new <- conv_w3$series
  
  expect_equal(res_lead.x$xpred_lead.new, conv_w$series)
  expect_equal(res_lead.x$xpred_targ.new, conv_w3$series)
})

test_that("output() returns a KFS object", {
  res <- make_gompertz_fit()
  kfs <- output(res)
  
  expect_true(inherits(kfs, "KFS"))
})

test_that("modelKFS() extracts SSModel from KFS", {
  res <- make_gompertz_fit()
  kfs <- output(res)
  mod <- modelKFS(kfs)
  
  expect_true(is.SSModel(mod))
})

test_that("seasonalComp() returns seasonal component info", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos)
  res <- estimate(model)
  
  sc <- seasonalComp(output(res))
  
  expect_true(sc == 3)
})

test_that("att() returns filtered state estimates", {
  res <- make_gompertz_fit()
  a <- att(output(res))
  
  expect_true(is.matrix(a) || is.array(a))
})

test_that("Ptt() returns filtered covariance matrices", {
  res <- make_gompertz_fit()
  P <- Ptt(output(res))
  
  expect_true(is.array(P))
})

test_that("get_V() returns smoothed covariance matrices", {
  res <- make_gompertz_fit()
  V <- get_V(output(res))
  
  expect_true(is.array(V))
})

test_that("matrixKFS() extracts model matrices", {
  res <- make_gompertz_fit()
  Z <- matrixKFS(output(res), "Z")
  
  expect_true(is.array(Z))
})

test_that("alphahat() returns smoothed state estimates", {
  res <- make_gompertz_fit()
  a_hat <- alphahat(output(res))
  
  expect_true(is.matrix(a_hat) || is.array(a_hat))
})

test_that("gety() extracts observation series from SSModel", {
  res <- make_gompertz_fit()
  y <- gety(modelKFS(output(res)))
  
  expect_true(is.numeric(y))
  expect_length(y, length(res$data))
})

test_that("gety.hat() extracts predictions from predict_all output", {
  res <- make_gompertz_fit()
  preds <- res$predict_all(n.ahead = 7)
  yhat <- gety.hat(preds)
  
  expect_true(is_idx_series(yhat))
  expect_equal(colnames(idx_values(yhat)), c("y.hat", "y.hat.upr", "y.hat.lwr"))
  expect_equal(length(yhat), 7)
  expect_true(is.numeric(idx_values(yhat)[, "y.hat"]))
  expect_true(is.numeric(idx_values(yhat)[, "y.hat.upr"]))
  expect_true(is.numeric(idx_values(yhat)[, "y.hat.lwr"]))
})

test_that("estimate() returns FilterResults object", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series,
    q = 0.005,
    end = end.pos
  )
  
  res <- estimate(model)
  
  expect_true(inherits(res, "FilterResults"))
})

test_that("SSModelDynamicGompertz print/summary/plot work", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-06"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = conv$calendar)
  
  expect_error(print(model), NA)
  expect_error(summary(model), NA)
  expect_error(plot(model), NA)
})

test_that("FilterResults print and summary work", {
  res <- make_gompertz_fit()
  
  expect_error(print(res), NA)
  expect_error(summary(res), NA)
})

test_that("SSModelLeadingIndicator print/summary/plot and FilterResultsLI print/summary work", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  start.pos <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  end.pos   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  
  out_eng <- SSModelLeadingIndicator(
    Y = conv$series,
    n.lag = 4,
    sea.period = 7,
    start = start.pos,
    end   = end.pos,
    calendar = conv$calendar
  )
  
  expect_error(print(out_eng), NA)
  expect_error(summary(out_eng), NA)
  expect_error(plot(out_eng), NA)
  
  res_eng <- estimate(out_eng)
  
  expect_error(print(res_eng), NA)
  expect_error(summary(res_eng), NA)
})