library(KFAS)
library(xts)
library(stats)
library(timetk)
library(zoo)

test_that("tsgc produces same LI output as KFAS", {
  set.seed(123)
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  tsgc_est <- tsgc_mod$estimate()
  
  y<-add_daily_ldl(ukitaly, LeadIndCol=1)
  
  y$newCases = stats::lag(y$newCases,n.lag)
  y$LDLcases = stats::lag(y$LDLcases,n.lag)
  y$cCases = stats::lag(y$cCases,n.lag)
  
  y[is.infinite(y)] <- NA
  
  y <- get_timeframe(na.omit(y),est.start)
  
  data_ldl <- get_timeframe(y, est.start, est.end)[,c("LDLcases","LDLhosp")]
  
  data_mat <- as.matrix(data_ldl)
  
  kfas_mod <- SSModel(data_mat ~ SSMtrend(degree = 2, 
                                                        Q = matrix(c(0,0,0,NA),2,2),
                                                        type = 'common') +
                        SSMtrend(degree = 1, Q = matrix(NA),index=1),
                      H = matrix(c(NA,0,0,NA),2,2))
  npar <- sum(is.na(kfas_mod$Q)) + sum(is.na(kfas_mod$H))
  kfas_fit <- fitSSM(kfas_mod, rep(0,npar))
  kfas_est <- KFS(kfas_fit$model)
  
  expect_equal(unname(as.matrix(tsgc_est$output$alphahat)), 
                      unname(as.matrix(kfas_est$alphahat)))
})

test_that("tsgc produces same LI output as KFAS with fixed q", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = 0.005,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  tsgc_est <- tsgc_mod$estimate()
  
  updatesn=function(pars, model, snr, order, index){
    if(any(is.na(model$Q))){
      Q <- as.matrix(model$Q[,,1])
      naQd  <- which(is.na(diag(Q)))
      naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))
      Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
      diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
      Q[naQd,naQd][naQnd] <- pars[length(naQd)+1:length(naQnd)]
      model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
    }
    if(!identical(model$H,'Omitted') && any(is.na(model$H))){
      H<-as.matrix(model$H[,,1])
      naHd  <- which(is.na(diag(H)))
      naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
      H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0
      diag(H)[naHd] <-
        exp(0.5 * pars[length(naQd)+length(naQnd)+1:length(naHd)])
      H[naHd,naHd][naHnd] <-
        pars[length(naQd)+length(naQnd)+length(naHd)+1:length(naHnd)]
      model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
      model$Q[order,order,1] <- snr*crossprod(H[index,index])
    }
    model
  }
  updateli = updatesn %>% purrr::partial(snr=0.005,order=2,index=2)
  
  y<-add_daily_ldl(ukitaly, LeadIndCol=1)
  
  y$newCases = stats::lag(y$newCases,n.lag)
  y$LDLcases = stats::lag(y$LDLcases,n.lag)
  y$cCases = stats::lag(y$cCases,n.lag)
  
  y[is.infinite(y)] <- NA
  
  y <- get_timeframe(na.omit(y),est.start)
  
  data_ldl <- get_timeframe(y, est.start, est.end)[,c("LDLcases","LDLhosp")]
  
  data_mat <- as.matrix(data_ldl)
  
  kfas_mod <- SSModel(data_mat ~ SSMtrend(degree = 2, 
                                                        Q = matrix(c(0,0,0,NA),2,2),
                                                        type = 'common') +
                        SSMtrend(degree = 1, Q = matrix(NA),index=1),
                      H = matrix(c(NA,0,0,NA),2,2))
  npar <- sum(is.na(kfas_mod$Q)) + sum(is.na(kfas_mod$H))
  kfas_fit <- fitSSM(kfas_mod, rep(0,npar), updatefn = updateli)
  kfas_est <- KFS(kfas_fit$model)
  
  tsgc_snr <- tsgc_est$output$model$Q[2,2,1]/tsgc_est$output$model$H[2,2,1]
  kfas_snr <- kfas_est$model$Q[2,2,1]/kfas_est$model$H[2,2,1]
  
  expect_equal(tsgc_snr, 0.005)
  expect_equal(kfas_snr, 0.005)
  expect_equal(unname(as.matrix(tsgc_est$output$alphahat)), 
               unname(as.matrix(kfas_est$alphahat)))
})

test_that("Summary method works", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  
  expect_no_error(expect_no_warning(tsgc_mod$summary()))
})

test_that("Print method works", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  
  expect_no_error(expect_no_warning(tsgc_mod$print()))
})

test_that("Plot method works", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  
  expect_no_error(suppressWarnings(tsgc_mod$plot()))
})

test_that("Leading indicator model + seasonal has correct number of elements", {
  eng <- tsgc::england[, 1:2]
  
  est.start.eng <- as.Date("2021-04-30")
  est.end.eng   <- as.Date("2021-07-24")
  sea <- 7
  
  mod <- SSModelLeadingIndicator$new(eng, n.lag = 5, start.date = est.start.eng, 
                                     end.date = est.end.eng, sea.period = sea)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat), 3 + 2*(sea-1))
})

test_that("LI model + xpred1 + seasonal has correct number of elements", {
  eng <- tsgc::england[, 1:2]
  
  est.start.eng <- as.Date("2021-04-30")
  est.end.eng   <- as.Date("2021-07-24")
  sea <- 7
  
  xp <- england_weather_2021[, 1:4]
  
  mod <- SSModelLeadingIndicator$new(eng, n.lag = 5, start.date = est.start.eng, 
                                     end.date = est.end.eng, sea.period = sea,
                                     xpred1 = xp)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat),3 + 2*(sea-1) + ncol(xp))
})

test_that("LI model + xpred2 has correct number of elements", {
  eng <- tsgc::england[, 1:2]
  
  est.start.eng <- as.Date("2021-04-30")
  est.end.eng   <- as.Date("2021-07-24")
  sea <- 0
  
  xp <- england_weather_2021[, 1:4]
  
  mod <- SSModelLeadingIndicator$new(eng, n.lag = 5, start.date = est.start.eng, 
                                     end.date = est.end.eng, sea.period = sea,
                                     xpred2 = xp)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat),3 + ncol(xp))
})

test_that("LI + xpred1 + xpred2 + seasonal has correct number of elements", {
  eng <- tsgc::england[, 1:2]
  
  est.start.eng <- as.Date("2021-04-30")
  est.end.eng   <- as.Date("2021-07-24")
  sea <- 7
  
  xp <- england_weather_2021[, 1:4]
  
  mod <- SSModelLeadingIndicator$new(eng, n.lag = 5, start.date = est.start.eng, 
                                     end.date = est.end.eng, sea.period = sea,
                                     xpred1 = xp, xpred2 = xp)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat),3 + 2*(sea-1) + 2*ncol(xp))
})

test_that("LI with quarterly data has correct number of components", {
  data(nintendo_sales, package = "tsgc")
  
  sea <- 4
  
  est.start.q2  <- zoo::as.yearqtr("2017 Q1")
  est.end.q2    <- zoo::as.yearqtr("2019 Q4")
  n.lag.q       <- zoo::as.yearqtr("2017 Q1") - zoo::as.yearqtr("2006 Q4")
  
  y_q <- nintendo_sales[, c("wii", "switch_all")]
  
  mod_switch <- tsgc::SSModelLeadingIndicator$new(
    Y = y_q, sea.period = sea, n.lag = n.lag.q,
    start.date = est.start.q2, end.date = est.end.q2
  )
  res <- mod_switch$estimate()
  
  expect_equal(ncol(res$output$alphahat),3 + 2*(sea-1))
})

test_that("LI with monthly data has correct number of components", {
  data(etrading_apps, package = "tsgc")
  
  sea <- 12
  
  est.start.m2 <- zoo::as.yearmon(2017.5)
  est.end.m2   <- zoo::as.yearmon(2021 + 1/12)
  n.lag.m      <- zoo::as.yearmon(2017.5) - zoo::as.yearmon(2017)
  
  y_m <- etrading_apps[, c("DEGIRO", "AvaTrade")]
  mod_500_lead <- tsgc::SSModelLeadingIndicator$new(
    Y = y_m, sea.period = 12, n.lag = n.lag.m,
    start.date = est.start.m2, end.date = est.end.m2
  )
  res <- mod_500_lead$estimate()
  
  expect_equal(ncol(res$output$alphahat),3 + 2*(sea-1))
})

test_that("LI with annual data has correct number of components", {
  data(nintendo_sales, package = "tsgc")
  est.start.y <- zoo::as.yearmon(2011)
  est.end.y   <- zoo::as.yearmon(2018)
  
  yearly_nintendo      <- nintendo_sales[4 * (1:19), c("wii", "3ds")]
  yearly_nintendo_xts  <- xts::xts(zoo::coredata(yearly_nintendo), 
                                   order.by = zoo::yearmon(2005:2023))
  
  est.start.m2 <- zoo::as.yearmon(2017.5)
  est.end.m2   <- zoo::as.yearmon(2021 + 1/12)
  n.lag.y <- zoo::as.yearmon(2011) - zoo::as.yearmon(2007)
  
  mod_lead_y <- tsgc::SSModelLeadingIndicator$new(
    Y = yearly_nintendo_xts, sea.period = 0, n.lag = n.lag.y,
    start.date = est.start.y, end.date = est.end.y, LeadIndCol = 1
  )
  res <- mod_lead_y$estimate()
  
  expect_equal(ncol(res$output$alphahat),3)
})