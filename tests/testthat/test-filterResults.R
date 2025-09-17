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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
})

test_that("predict_level computes predictions correctly - seasonal + xpred + *NO* AR1, sea.on = TRUE", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, xpred = gauteng_weather, sea.period = 7,
    start.date = est.start, end.date = est.end
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
                       + SSMregression(~new_weather), 
                       H = matrix(Ht))
  
  delta_pred <- predict(res$output$model, newdata = new_model, 
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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
})

test_that("predict_level computes predictions correctly - seasonal + xpred + *NO* AR1, sea.on = FALSE", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start <- as.Date("2021-02-01")
  est.end  <- as.Date("2021-04-19")
  nf <- 7
  
  model <- SSModelDynamicGompertz$new(
    Y = gauteng$cum_cases, xpred = gauteng_weather, sea.period = 7,
    start.date = est.start, end.date = est.end
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
                       + SSMregression(~new_weather), 
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
  expect_equal(unname(as.vector(forc_tsgc$lower)), forc_lwr, tolerance = 1)
  expect_equal(unname(as.vector(forc_tsgc$upper)), forc_upr, tolerance = 1)
})

test_that("Test predict_all() gives same results as KFAS", {
  data('gauteng')
  q <- 0.005
  n.ahead <- 14
  y <- gauteng[1:50]
  estimation.date.start <- index(y)[1]

  y.new <- gauteng[1:64]
  y.new[51:64] <- NA
  model <- tsgc::SSModelDynamicGompertz$new(Y = y, q = q)
  res <- model$estimate()
  model <- tsgc::SSModelDynamicGompertz$new(Y=y.new, q = q)
  res.kfas <- model$estimate()
  filtered.out <- res$predict_all(n.ahead = 14, return.all = TRUE,
                                  sea.on = TRUE)

  expect_true(all(res.kfas$output$att == filtered.out$a.t.t))
  expect_equal(res.kfas$output$Ptt, filtered.out$P.t.t)

  # A solution - Extract model - extend y, change n ensure int then do the
  # KFS again.
  new.model <- res$output$model
  new.model$y <- rbind(
    new.model$y,
    matrix(NA, ncol = ncol(new.model$y), nrow = n.ahead) %>% as.ts()
  )
  attr(new.model, 'n') <- 64 %>% as.integer()
  model_output <- KFS(new.model)
  expect_equal(model_output$Ptt, filtered.out$P.t.t, tolerance = 1e-10)

  Zt <- drop(res.kfas$output$model$Z)
  Tt <- drop(res.kfas$output$model$T)
  y.t.t <- drop(Zt %*% t(res.kfas$output$att))
  expect_equal(y.t.t[51:64] %>% as.numeric(), filtered.out$y[51:64] %>%
    as.numeric())

  # Can use predict on non-extended y.
  y.hat.kfas <- predict(
    res$output$model, interval = c('prediction'), n.ahead = n.ahead,
    level = 0.68, states = c('all')
  )
  dates <- seq(tail(res.kfas$index,1) + 1, by = 'day', length.out = n.ahead)
  y.hat.kfas <- xts(y.hat.kfas[,1], order.by = dates)
  expect_equal(y.hat.kfas %>% as.numeric(), filtered.out$y[51:64] %>%
    as.numeric(), tolerance = 1e-10)

})



