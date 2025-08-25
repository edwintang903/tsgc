library(KFAS)

test_that("tsgc gives same output as KFAS", {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  dat.ldl <- na.omit(df2ldl(dat))
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                       q = NULL,
                                       sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  final.states.tsgc <- tsgc.est$output$alphahat
  
  kfas.model <- SSModel(as.matrix(dat.ldl) ~
                                SSMtrend(degree = 2, 
                                               Q = list(matrix(0), matrix(NA))), 
                              H = matrix(NA))
  kfas.fit <- fitSSM(kfas.model,inits=c(0,0))
  kfas.out <- KFS(kfas.fit$model)
  final.states.kfas <- kfas.out$alphahat
  
  expect_equal(final.states.tsgc,final.states.kfas,tolerance=1e-4)
})

test_that({"tsgc enforces signal-to-noise restrictions correctly"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  
  q.out <- as.vector(tsgc.est$output$model$Q[2,2,1]/tsgc.est$output$model$H)
  
  expect_equal(q.out,q.choose)
})

test_that({"tsgc summary method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.model$summary()
})

test_that({"tsgc print method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.model$print()
})

test_that({"tsgc plot method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  suppressWarnings(tsgc.model$plot())
})

test_that({"tsgc enforces signal-to-noise restrictions correctly"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  
  q.out <- as.vector(tsgc.est$output$model$Q[2,2,1]/tsgc.est$output$model$H)
  
  expect_equal(q.out,q.choose)
})


test_that({"Model with seasonal has expected number of seasonal components"},{
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  sea.choose <- 7
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = 0.005,
                                           sea.period = sea.choose)
  tsgc.est <- tsgc.model$estimate()
  sea.est <- length(tsgc.est$output$model["a1","seasonal"])
  
  expect_equal(sea.est,2*floor(sea.choose/2))
})

#test_that({"Model with xpred has expected number of slope coefficients"},{
#  data(gauteng, package = "tsgc")
#  cumulative_cases <- gauteng[, 1]
#  data(gauteng_weather_2021, package = "tsgc")
#  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
#  
#  est.start.1 <- as.Date("2021-02-01")
#  est.end.1   <- as.Date("2021-04-19")
#  
#  model_weather <- tsgc::SSModelDynamicGompertz$new(
#    Y = cumulative_cases,
#    xpred = gauteng_weather,
#    start.date = est.start.1,
#    end.date   = est.end.1,
#    sea.period = 0
#  )
#
#  est_weather <- model_weather$estimate()
#
#})

#test_that("Reinitialised model uses prior information correctly", {
#  data(gauteng, package = "tsgc")
#  cumulative_cases <- gauteng[, 1]
#  est.start.1 <- as.Date("2021-02-01")
#  est.end.2 <- as.Date("2021-06-25")
#  est.end.1   <- as.Date("2021-04-19")
#  reinit.date <- as.Date("2021-04-21")
#  
#  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
#    Y = cumulative_cases, q = q.default,
#    start.date = est.start.1, end.date = est.end.2
#  )
#  
#  model_reinit <- tsgc::SSModelDynamicGompertz$new(
#    Y = cumulative_cases, q = q.default,
#    start.date = est.start.1, end.date = est.end.2,
#    reinit.date = reinit.date
#  )
#})

test_that("tsgc gives same output as KFAS for fixed q", {
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
  
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  dat.ldl <- na.omit(df2ldl(dat))
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = 0.005,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  final.states.tsgc <- tsgc.est$output$alphahat
  
  kfas.model <- SSModel(as.matrix(dat.ldl) ~
                          SSMtrend(degree = 2, 
                                   Q = list(matrix(0), matrix(NA))), 
                        H = matrix(NA))
  npar <- sum(is.na(kfas.model$Q)) + sum(is.na(kfas.model$H))
  update <- purrr::partial(updatesn, snr=0.005, order = 2, index = 1)
  kfas.fit <- fitSSM(kfas.model,inits=rep(0,npar), updatefn = update)
  kfas.out <- KFS(kfas.fit$model)
  final.states.kfas <- kfas.out$alphahat
  
  expect_equal(final.states.tsgc,final.states.kfas,tolerance=1e-4)
})
