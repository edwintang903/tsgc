library(tsgc)
library(zoo)
library(KFAS)

make_gompertz_fit <- function(end.date = as.Date("2020-07-20")) {
  model <- SSModelDynamicGompertz$new(
    Y = gauteng,
    q = 0.005,
    end.date = end.date
  )
  estimate(model)
}

test_that("supply_xpred.new supplies future xpred values to FilterResults object", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather<-gauteng_weather_2021[,c(1,3)]
  
  # Set up model and estimate it
  model_weather <- SSModelDynamicGompertz$new(Y = gauteng, xpred=gauteng_weather,
                                              start.date=as.Date("2021-02-01"), 
                                              end.date=as.Date("2021-04-19"))
  res_weather <- estimate(model_weather)
  res_weather$xpred.new
  
  # Feed future weather data into the results object. Subsetting of gauteng_weather 
  #is done inside the function.
  supply_xpred.new(res_weather,gauteng_weather)
  expect_equal(res_weather$xpred.new, gauteng_weather)
})

test_that("supply_xpred.new supplies future xpred values to FilterResultsLI object", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  
  # Set up model and estimate it
  mod<-SSModelLeadingIndicator$new(england[,1:2], n.lag=4, 
                                   xpred_lead=england_weather_2021[,1:2], 
                                   xpred_targ=england_weather_2021[,3], 
                                   start.date = as.Date("2021-04-30"), 
                                   end.date = as.Date("2021-07-24"))
  res_lead.x<-estimate(mod)
  
  supply_xpred.new(res_lead.x,england_weather_2021[,1:2],idx='lead')
  supply_xpred.new(res_lead.x,england_weather_2021[,3],idx='targ')
  
  expect_equal(res_lead.x$xpred_lead.new, england_weather_2021[,1:2])
  expect_equal(res_lead.x$xpred_targ.new, england_weather_2021[,3])
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
  idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
  model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
  res <- estimate(model)
  
  sc <- seasonalComp(output(res))
  
  expect_true(sc==3)
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
  expect_length(y, length(res$data_xts))
})


test_that("gety.hat() extracts predictions from predict_all output", {
  res <- make_gompertz_fit()
  preds <- res$predict_all(n.ahead = 7)
  yhat <- gety.hat(preds)
  
  # Structure
  expect_true(is.matrix(yhat))
  
  # Required columns
  expect_true(all(c("y.hat", "y.hat.upr", "y.hat.lwr") %in% colnames(yhat)))
  
  # Dimensions
  expect_equal(nrow(yhat), 7)
  expect_equal(ncol(yhat), 3)
  
  # Column types
  expect_true(is.numeric(yhat[, "y.hat"]))
  expect_true(is.numeric(yhat[, "y.hat.upr"]))
  expect_true(is.numeric(yhat[, "y.hat.lwr"]))
})

test_that("estimate() returns FilterResults object", {
  model <- SSModelDynamicGompertz$new(
    Y = gauteng,
    q = 0.005,
    end.date = as.Date("2020-07-20")
  )
  
  res <- estimate(model)
  
  expect_true(inherits(res, "FilterResults"))
})

test_that("SSModelDynamicGompertz print/summary/plot work", {
  idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
  model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
  
  expect_error(print(model), NA)
  expect_error(summary(model), NA)
  expect_error(plot(model), NA)
})

test_that("FilterResults print and summary work", {
  res <- make_gompertz_fit()
  
  expect_error(print(res),NA)
  expect_error(summary(res),NA)
})

test_that("SSModelLeadingIndicator print/summary/plot 
          and FilterResultsLI print/summary work", {
  out_eng <- SSModelLeadingIndicator(
    Y = england[, 1:2],
    n.lag = 4,
    sea.period = 7,
    start.date = as.Date("2021-04-30"),
    end.date   = as.Date("2021-07-24")
  )
  
  expect_error(print(out_eng),NA)
  expect_error(summary(out_eng),NA)
  expect_error(plot(out_eng),NA)
  
  res_eng <- estimate(out_eng)
  
  expect_error(print(res_eng),NA)
  expect_error(summary(res_eng),NA)
})
