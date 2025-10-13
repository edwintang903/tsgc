library(KFAS)

test_that("predict_level computes predictions correctly - no seasonal", {
  data(gauteng, package = "tsgc")
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, q = 0.005, sea.period = 0,
    start.date = est.start, end.date = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf, 
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions of cumulated variable correctly - no seasonal", {
  data(gauteng, package = "tsgc")
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, q = 0.005, sea.period = 0,
    start.date = est.start, end.date = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf, 
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  mult <- cumprod(1+exp(delta_fit))
  forc <- rep(YT,nf)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  mult_lwr <- cumprod(1+exp(delta_lwr))
  forc_lwr <- rep(YT,nf)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, return.diff = FALSE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
})

test_that("predict_level computes predictions correctly - seasonal, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, q = 0.005, sea.period = 7,
    start.date = est.start, end.date = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf, 
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = TRUE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal but sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, q = 0.005, sea.period = 7,
    start.date = est.start, end.date = est.end
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf, 
                        interval = c("confidence"), level = 0.68,
                        states = c("level"))
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + AR1, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, q = 0.005, sea.period = 7,
    start.date = est.start, end.date = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf, 
                        interval = c("confidence"), level = 0.68)
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = TRUE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + AR1, sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, q = 0.005, sea.period = 7,
    start.date = est.start, end.date = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  delta_pred <- predict(res$output$model, n.ahead = nf, 
                        interval = c("confidence"), level = 0.68,
                        states = c("level", "custom"))
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + xpred + AR1, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, xpred = gauteng_weather, sea.period = 7,
    start.date = est.start, end.date = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  supply_xpred.new(res, gauteng_weather)
  
  f.start <- est.end + 1
  f.end <- est.end + nf
  
  new_weather <- get_timeframe(gauteng_weather, f.start, f.end)
  
  Qt.slope <- res$output$model$Q[2,2,1]
  Qt.seas <- res$output$model$Q[3,3,1]
  Qt.ar1 <- res$output$model$Q[9,9,1]
  Ht <- res$output$model$H[1,1,1]
  
  new_model <- SSModel(formula = matrix(rep(NA,nf), ncol = 1) ~ 
                         SSMtrend(degree = 2, Q = list(matrix(0), 
                                                       matrix(Qt.slope))) 
                       + SSMseasonal(period = 7, Q = Qt.seas,
                                     sea.type = "trigonometric") 
                       + SSMregression(~new_weather)
                       + SSMcustom(Z=1,T=1,R=1,Q=Qt.ar1,
                                   state_names="ar1"), 
                       H = matrix(Ht))
  
  delta_pred <- predict(res$output$model, newdata = new_model, 
                        interval = c("confidence"), level = 0.68, states = c("all"))
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = TRUE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level computes predictions correctly - seasonal + xpred + AR1, sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, xpred = gauteng_weather, sea.period = 7,
    start.date = est.start, end.date = est.end, ar1 = TRUE
  )
  res <- estimate(model)
  
  supply_xpred.new(res, gauteng_weather)
  
  f.start <- est.end + 1
  f.end <- est.end + nf
  
  new_weather <- get_timeframe(gauteng_weather, f.start, f.end)
  
  Qt.slope <- res$output$model$Q[2,2,1]
  Qt.seas <- res$output$model$Q[3,3,1]
  Ht <- res$output$model$H[1,1,1]
  
  new_model <- SSModel(formula = matrix(rep(NA,nf), ncol = 1) ~ 
                         SSMtrend(degree = 2, Q = list(matrix(0), 
                                                       matrix(Qt.slope))) 
                       + SSMseasonal(period = 7, Q = Qt.seas,
                                     sea.type = "trigonometric") 
                       + SSMregression(~new_weather)
                       + SSMcustom(Z=1,T=1,R=1,Q=Qt.ar1,
                                   state_names="ar1"), 
                       H = matrix(Ht))
  
  delta_pred <- predict(res$output$model, newdata = new_model, 
                        interval = c("confidence"), level = 0.68,
                        states = c("level", "custom", "regression"))
  
  delta_fit <- as.vector(delta_pred[,"fit"])
  YT <- tail(model$Y,1)
  cp <- cumprod(1+exp(delta_fit))
  mult <- c(1,cp[1:(nf-1)])
  forc <- rep(YT,nf)*exp(delta_fit)*mult
  
  delta_lwr <- as.vector(delta_pred[,"lwr"])
  cp_lwr <- cumprod(1+exp(delta_lwr))
  mult_lwr <- c(1,cp_lwr[1:(nf-1)])
  forc_lwr <- rep(YT,nf)*exp(delta_lwr)*mult_lwr
  
  delta_upr <- as.vector(delta_pred[,"upr"])
  cp_upr <- cumprod(1+exp(delta_upr))
  mult_upr <- c(1,cp_upr[1:(nf-1)])
  forc_upr <- rep(YT,nf)*exp(delta_upr)*mult_upr
  
  forc_tsgc <- res$predict_level(n.ahead = nf, sea.on = FALSE)
  
  expect_equal(unname(as.vector(forc_tsgc$fit)), forc)
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 0.005)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 0.005)
})

test_that("predict_level works - quarterly data", {
  data(nintendo_sales, package = "tsgc")
  wii <- nintendo_sales[, 1]
  
  est.start.q <- zoo::as.yearqtr("2006 Q4")
  est.end.q   <- zoo::as.yearqtr("2010 Q3")
  
  nf <- 4
  
  mod_wii <- tsgc::SSModelDynamicGompertz$new(
    Y = wii, sea.period = 4, start.date = est.start.q, end.date = est.end.q
  )
  res_wii <- mod_wii$estimate()
  
  forecasts <- res_wii$predict_level(n.ahead = nf)
  
  expect_equal(length(forecasts$fit), nf)
})

test_that("predict_level works - monthly data", {
  data(etrading_apps, package = "tsgc")
  Plus500 <- etrading_apps[, 1]
  
  est.start.m <- zoo::as.yearmon(2016)
  est.end.m   <- zoo::as.yearmon(2021)
  
  nf <- 12
  
  mod_500 <- tsgc::SSModelDynamicGompertz$new(
    Y = Plus500, sea.period = 12, start.date = est.start.m, end.date = est.end.m
  )
  res_500 <- mod_500$estimate()
  
  forecasts <- res_500$predict_level(n.ahead = nf)
  
  expect_equal(length(forecasts$fit), nf)
})

test_that("predict_level works - annual data", {
  data(nintendo_sales, package = "tsgc")
  
  est.start.y <- zoo::as.yearmon(2011)
  est.end.y   <- zoo::as.yearmon(2018)
  
  yearly_nintendo      <- nintendo_sales[4 * (1:19), c("wii", "3ds")]
  threeds_xts          <- xts::xts(zoo::coredata(yearly_nintendo[, "3ds"]), order.by = zoo::yearmon(2005:2023))
  yearly_nintendo_xts  <- xts::xts(zoo::coredata(yearly_nintendo), order.by = zoo::yearmon(2005:2023))
  
  mod_3ds <- tsgc::SSModelDynamicGompertz$new(
    Y = threeds_xts, sea.period = 0, start.date = est.start.y, end.date = est.end.y
  )
  
  res_3ds <- estimate(mod_3ds)
  
  forecasts <- res_3ds$predict_level(n.ahead = 1)
  
  expect_equal(length(forecasts$fit),1)
})


