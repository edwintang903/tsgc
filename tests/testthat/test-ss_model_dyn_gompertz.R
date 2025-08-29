library(KFAS)
library(xts)

test_that("tsgc gives same output as KFAS", {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  dat.ldl <- df2ldl(dat)
  
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
  dat.ldl <- df2ldl(dat)
  
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

test_that({"tsgc summary method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  expect_no_error(expect_no_warning(tsgc.model$summary()))
})

test_that({"tsgc print method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  expect_no_error(expect_no_warning(tsgc.model$print()))
})

test_that({"tsgc plot method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  expect_no_error(suppressWarnings(tsgc.model$plot()))
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
  
  expect_equal(sea.est,(sea.choose-1))
})

test_that({"Model with xpred has expected number of slope coefficients"},{
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start.1 <- as.Date("2021-02-01")
  est.end.1   <- as.Date("2021-04-19")
  
  model_weather <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases,
    xpred = gauteng_weather,
    start.date = est.start.1,
    end.date   = est.end.1,
    sea.period = 0
  )

  est_weather <- model_weather$estimate()
  
  expect_equal(length(est_weather$output$model["a1","regression"]),
               ncol(gauteng_weather))
})

test_that({"Model with xpred and seasonal has expected number of 
  seasonal and slope coefficients"},{
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start.1 <- as.Date("2021-02-01")
  est.end.1   <- as.Date("2021-04-19")
  
  sea.choose = 7
  
  model_weather <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases,
    xpred = gauteng_weather,
    start.date = est.start.1,
    end.date   = est.end.1,
    sea.period = sea.choose
  )
  
  est_weather <- model_weather$estimate()
  
  expect_equal(length(est_weather$output$model["a1","regression"]),
               ncol(gauteng_weather))
  expect_equal(length(est_weather$output$model["a1","seasonal"]),
               (sea.choose-1))
})

test_that("Reinitialised model uses prior information correctly", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  est.start.1 <- as.Date("2021-02-01")
  est.end.2 <- as.Date("2021-06-25")
  reinit.date <- as.Date("2021-04-21")
  q.default <- NULL
  sea.default <- 0
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default,
    start.date = est.start.1, end.date = reinit.date
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default,
    start.date = est.start.1, end.date = est.end.2,
    reinit.date = reinit.date
  )
  
  i.reinit <- which(index(model_rei_base$Y)==reinit.date)
  lYy <- as.numeric(
    log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit,]
  Tt <- est_base$output$model$T[,,1]
  Ptt <- est_base$output$Ptt[,,i.reinit]
  Rt <- est_base$output$model$R[,,1]
  Qt <- est_base$output$model$Q[,,1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  level.idx <- which(names(att)=="level")
  slope.idx <- which(names(att)=="slope")
  
  adj <- rep(0, length(att))
  adj[level.idx] <- lYy
  
  a1.base <- Tt %*% att + adj
  a1.base[slope.idx] = 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base, a1.reinit)
  expect_equal(P1.base, P1.reinit)
})

test_that("Reinitialised model with seasonal uses prior information correctly", 
  {
    data(gauteng, package = "tsgc")
    cumulative_cases <- gauteng[, 1]
    est.start.1 <- as.Date("2021-02-01")
    est.end.2 <- as.Date("2021-06-25")
    reinit.date <- as.Date("2021-04-21")
    q.default <- NULL
    sea.default <- 7
    
    model_rei_base <- tsgc::SSModelDynamicGompertz$new(
      Y = cumulative_cases, q = q.default, sea.period = sea.default,
      start.date = est.start.1, end.date = reinit.date
    )
    
    model_reinit <- tsgc::SSModelDynamicGompertz$new(
      Y = cumulative_cases, q = q.default, sea.period = sea.default,
      start.date = est.start.1, end.date = est.end.2,
      reinit.date = reinit.date
    )
    
    i.reinit <- which(index(model_rei_base$Y)==reinit.date)
    lYy <- as.numeric(
      log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
    
    est_base <- model_rei_base$estimate()
    est_reinit <- model_reinit$estimate()
    
    att <- est_base$output$att[i.reinit,]
    Tt <- est_base$output$model$T[,,1]
    Ptt <- est_base$output$Ptt[,,i.reinit]
    Rt <- est_base$output$model$R[,,1]
    Qt <- est_base$output$model$Q[,,1]
    
    P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
    
    a1.base <- Tt %*% att
    level.idx <- which(rownames(a1.base)=="level")
    slope.idx <- which(rownames(a1.base)=="slope")
    sea.rows <- grep("sea_trig",rownames(a1.base))
    adj <- rep(0, length(a1.base))
    adj[level.idx] <- lYy
    a1.base <- a1.base + adj
    a1.base[slope.idx] <- 0
    a1.reinit <- est_reinit$output$model$a1
    
    P1.reinit <- est_reinit$output$model$P1
    
    expect_equal(a1.base, a1.reinit)
    expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
                 P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
    expect_equal(P1.base[sea.rows,sea.rows],P1.reinit[sea.rows,sea.rows])
  }
)

test_that("Reinitialised model with xpred uses prior information correctly", {
      data(gauteng, package = "tsgc")
      cumulative_cases <- gauteng[, 1]
      est.start.1 <- as.Date("2021-02-01")
      est.end.2 <- as.Date("2021-06-25")
      reinit.date <- as.Date("2021-04-21")
      q.default <- NULL
      sea.default <- 0
      
      data(gauteng_weather_2021, package = "tsgc")
      gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
      
      model_rei_base <- tsgc::SSModelDynamicGompertz$new(
        Y = cumulative_cases, q = q.default, sea.period = sea.default,
        xpred = gauteng_weather, start.date = est.start.1, end.date = reinit.date
      )
      
      model_reinit <- tsgc::SSModelDynamicGompertz$new(
        Y = cumulative_cases, q = q.default, sea.period = sea.default,
        xpred = gauteng_weather, start.date = est.start.1, end.date = est.end.2,
        reinit.date = reinit.date
      )
      
      i.reinit <- which(index(model_rei_base$Y)==reinit.date)
      lYy <- as.numeric(
        log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
      
      est_base <- model_rei_base$estimate()
      est_reinit <- model_reinit$estimate()
      
      att <- est_base$output$att[i.reinit,]
      Tt <- est_base$output$model$T[,,1]
      Ptt <- est_base$output$Ptt[,,i.reinit]
      Rt <- est_base$output$model$R[,,1]
      Qt <- est_base$output$model$Q[,,1]
      
      P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
      
      a1.base <- Tt %*% att
      level.idx <- which(rownames(a1.base)=="level")
      slope.idx <- which(rownames(a1.base)=="slope")
      xpred.rows <- grep("xpred",rownames(a1.base))
      adj <- rep(0, length(a1.base))
      adj[level.idx] <- lYy
      a1.base <- a1.base + adj
      a1.base[c(slope.idx,xpred.rows)] <- 0
      a1.reinit <- est_reinit$output$model$a1
      
      P1.base[xpred.rows,xpred.rows] <- 0
      P1.reinit <- est_reinit$output$model$P1
      
      expect_equal(a1.base, a1.reinit)
      expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
                   P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
      expect_equal(P1.base[xpred.rows,xpred.rows],
                   P1.reinit[xpred.rows,xpred.rows])
    }
)

test_that("Reinitialised model with xpred and seasonal uses prior information correctly", {
    data(gauteng, package = "tsgc")
    cumulative_cases <- gauteng[, 1]
    est.start.1 <- as.Date("2021-02-01")
    est.end.2 <- as.Date("2021-06-25")
    reinit.date <- as.Date("2021-04-21")
    q.default <- NULL
    sea.default <- 7
    
    data(gauteng_weather_2021, package = "tsgc")
    gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
    
    model_rei_base <- tsgc::SSModelDynamicGompertz$new(
      Y = cumulative_cases, q = q.default, sea.period = sea.default,
      xpred = gauteng_weather, start.date = est.start.1, end.date = reinit.date
    )
    
    model_reinit <- tsgc::SSModelDynamicGompertz$new(
      Y = cumulative_cases, q = q.default, sea.period = sea.default,
      xpred = gauteng_weather, start.date = est.start.1, end.date = est.end.2,
      reinit.date = reinit.date
    )
    
    i.reinit <- which(index(model_rei_base$Y)==reinit.date)
    lYy <- as.numeric(
      log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
    
    est_base <- model_rei_base$estimate()
    est_reinit <- model_reinit$estimate()
    
    att <- est_base$output$att[i.reinit,]
    Tt <- est_base$output$model$T[,,1]
    Ptt <- est_base$output$Ptt[,,i.reinit]
    Rt <- est_base$output$model$R[,,1]
    Qt <- est_base$output$model$Q[,,1]
    
    P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
    
    a1.base <- Tt %*% att
    level.idx <- which(rownames(a1.base)=="level")
    slope.idx <- which(rownames(a1.base)=="slope")
    sea.rows <- grep("sea_trig",rownames(a1.base))
    xpred.rows <- grep("xpred",rownames(a1.base))
    adj <- rep(0, length(a1.base))
    adj[level.idx] <- lYy
    a1.base <- a1.base + adj
    a1.base[c(slope.idx,xpred.rows)] <- 0
    a1.reinit <- est_reinit$output$model$a1
    
    P1.base[xpred.rows,xpred.rows] <- 0
    P1.reinit <- est_reinit$output$model$P1
    
    expect_equal(a1.base, a1.reinit)
    expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
                 P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
    expect_equal(P1.base[sea.rows,sea.rows],P1.reinit[sea.rows,sea.rows])
    expect_equal(P1.base[xpred.rows,xpred.rows],P1.reinit[xpred.rows,xpred.rows])
  }
)

test_that("Vanilla AR(1) Gompertz curve gives same output as KFAS", {
  data('england', package = 'tsgc')
  daily_cases <- diff(england$cum_cases)[2:31]
  cumulative_cases <- england$cum_cases[2:31] #cumsum(daily_cases)
  
  # Create ldl series
  ldl_cases <- log(daily_cases/lag(cumulative_cases))
  
  # KFAS model
  mod <- SSModel(as.matrix(ldl_cases) ~
                   SSMtrend(degree = 2, Q = list(matrix(0), matrix(NA))) +
                   SSMcustom(Z = 1, T = 1, R = 1, Q = matrix(NA), state_names = "ar1"),
                 H = matrix(NA))
  
  npar <- sum(is.na(mod$Q)) + sum(is.na(mod$H))
  
  update_ar1 <- function(pars, model) {
    
    H <- as.matrix(model$H[,,1])
    H[1,1] <- exp(0.5 * pars[1])
    model$H[,,1] <- crossprod(H)
    
    Q <- as.matrix(model$Q[,,1])
    Q[2,2] <- exp(0.5 * pars[2])      
    Q[3,3] <- exp(0.5 * pars[3])      
    Q[2:3,2:3] <- crossprod(Q[2:3,2:3])
    model$Q[,,1] <- Q
    
    Tmat <- model$T[,,1]
    Tmat[3,3] <- pars[4]
    model$T[,,1] <- Tmat
    
    return(model)
  }
  
  fit <- fitSSM(mod, inits = c(0,0,0,1), updatefn = update_ar1)
  out <- KFS(fit$model)
  
  # tsgc model
  
  xts_cases = xts::xts(cumulative_cases,order.by=(as.Date('2021-01-01')+1:30))
  tsgc_mod <- SSModelDynamicGompertz$new(Y = xts_cases, ar1 = TRUE, sea.period = 0)
  tsgc_est <- tsgc_mod$estimate()
  
  expect_equal(out$alphahat[30,], tsgc_est$output$alphahat[30,],  tolerance = 1e-4)
})

test_that("Vanilla model with AR(1) has correct number of components", {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  dat.ldl <- na.omit(df2ldl(dat))
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = NULL,
                                           sea.period = 0,
                                           ar1 = TRUE)
  tsgc.est <- tsgc.model$estimate()
  final.states.tsgc <- tsgc.est$output$alphahat
  

  expect_equal(length(tsgc.est$output$alphahat[100,]),3)
})

test_that({"Model with seasonal + AR1 has expected number of components"},{
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  sea.choose <- 7
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = 0.005,
                                           sea.period = sea.choose,
                                           ar1 = TRUE)
  tsgc.est <- tsgc.model$estimate()

  expect_equal(length(tsgc.est$output$alphahat[100,]),(sea.choose-1)+3)
})

test_that({"Model with xpred and AR(1) has expected number of components"},{
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  est.start.1 <- as.Date("2021-02-01")
  est.end.1   <- as.Date("2021-04-19")
  
  model_weather <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases,
    xpred = gauteng_weather,
    start.date = est.start.1,
    end.date   = est.end.1,
    sea.period = 0,
    ar1 = TRUE
  )
  
  est_weather <- model_weather$estimate()
  
  expect_equal(length(est_weather$output$alphahat[10,]),
               ncol(gauteng_weather) + 3)
})

test_that({"Model with xpred, seasonal and AR(1) has expected number of components"},{
    data(gauteng, package = "tsgc")
    cumulative_cases <- gauteng[, 1]
    data(gauteng_weather_2021, package = "tsgc")
    gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
    
    est.start.1 <- as.Date("2021-02-01")
    est.end.1   <- as.Date("2021-04-19")
    
    sea.choose = 7
    
    model_weather <- tsgc::SSModelDynamicGompertz$new(
      Y = cumulative_cases,
      xpred = gauteng_weather,
      start.date = est.start.1,
      end.date   = est.end.1,
      sea.period = sea.choose,
      ar1 = TRUE
    )
    
    est_weather <- model_weather$estimate()
    
    expect_equal(length(est_weather$output$alphahat[10,]),
                 ncol(gauteng_weather) + (sea.choose-1) + 3)
  })

test_that("Reinitialised model with AR1 uses prior information correctly", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  est.start.1 <- as.Date("2021-02-01")
  est.end.2 <- as.Date("2021-06-25")
  reinit.date <- as.Date("2021-04-21")
  q.default <- NULL
  sea.default <- 0
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default,
    start.date = est.start.1, end.date = reinit.date, ar1 = TRUE
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default,
    start.date = est.start.1, end.date = est.end.2, ar1 = TRUE,
    reinit.date = reinit.date
  )
  
  i.reinit <- which(index(model_rei_base$Y)==reinit.date)
  lYy <- as.numeric(
    log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit,]
  Tt <- est_base$output$model$T[,,1]
  Ptt <- est_base$output$Ptt[,,i.reinit]
  Rt <- est_base$output$model$R[,,1]
  Qt <- est_base$output$model$Q[,,1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  level.idx <- which(names(att)=="level")
  slope.idx <- which(names(att)=="slope")
  ar1.idx <- which(names(att)=="ar1")
  
  adj <- rep(0, length(att))
  adj[level.idx] <- lYy
  
  a1.base <- Tt %*% att + adj
  a1.base[c(slope.idx,ar1.idx)] = 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base, a1.reinit)
  expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
               P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
})

test_that("Reinitialised model with seasonal and AR1 uses prior information correctly", 
          {
            data(gauteng, package = "tsgc")
            cumulative_cases <- gauteng[, 1]
            est.start.1 <- as.Date("2021-02-01")
            est.end.2 <- as.Date("2021-06-25")
            reinit.date <- as.Date("2021-04-21")
            q.default <- NULL
            sea.default <- 7
            
            model_rei_base <- tsgc::SSModelDynamicGompertz$new(
              Y = cumulative_cases, q = q.default, sea.period = sea.default,
              start.date = est.start.1, end.date = reinit.date, ar1 = TRUE
            )
            
            model_reinit <- tsgc::SSModelDynamicGompertz$new(
              Y = cumulative_cases, q = q.default, sea.period = sea.default,
              start.date = est.start.1, end.date = est.end.2, ar1 = TRUE,
              reinit.date = reinit.date
            )
            
            i.reinit <- which(index(model_rei_base$Y)==reinit.date)
            lYy <- as.numeric(
              log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
            
            est_base <- model_rei_base$estimate()
            est_reinit <- model_reinit$estimate()
            
            att <- est_base$output$att[i.reinit,]
            Tt <- est_base$output$model$T[,,1]
            Ptt <- est_base$output$Ptt[,,i.reinit]
            Rt <- est_base$output$model$R[,,1]
            Qt <- est_base$output$model$Q[,,1]
            
            P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
            
            a1.base <- Tt %*% att
            level.idx <- which(rownames(a1.base)=="level")
            slope.idx <- which(rownames(a1.base)=="slope")
            sea.rows <- grep("sea_trig",rownames(a1.base))
            choose <- c(level.idx,slope.idx,sea.rows)
            adj <- rep(0, length(a1.base))
            adj[level.idx] <- lYy
            a1.base <- a1.base + adj
            a1.base[slope.idx] <- 0
            a1.reinit <- est_reinit$output$model$a1
            
            P1.reinit <- est_reinit$output$model$P1
            
            expect_equal(a1.base[choose], a1.reinit[choose])
            expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
                         P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
            expect_equal(P1.base[sea.rows,sea.rows],P1.reinit[sea.rows,sea.rows])
          }
)

test_that("Reinitialised model with xpred and AR1 uses prior information correctly", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  est.start.1 <- as.Date("2021-02-01")
  est.end.2 <- as.Date("2021-06-25")
  reinit.date <- as.Date("2021-04-21")
  q.default <- NULL
  sea.default <- 0
  
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default, ar1 = TRUE,
    xpred = gauteng_weather, start.date = est.start.1, end.date = reinit.date
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default, ar1 = TRUE,
    xpred = gauteng_weather, start.date = est.start.1, end.date = est.end.2,
    reinit.date = reinit.date
  )
  
  i.reinit <- which(index(model_rei_base$Y)==reinit.date)
  lYy <- as.numeric(
    log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit,]
  Tt <- est_base$output$model$T[,,1]
  Ptt <- est_base$output$Ptt[,,i.reinit]
  Rt <- est_base$output$model$R[,,1]
  Qt <- est_base$output$model$Q[,,1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base)=="level")
  slope.idx <- which(rownames(a1.base)=="slope")
  xpred.rows <- grep("xpred",rownames(a1.base))
  choose <- c(level.idx, slope.idx, xpred.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[c(slope.idx,xpred.rows)] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.base[xpred.rows,xpred.rows] <- 0
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
               P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
  expect_equal(P1.base[xpred.rows,xpred.rows],
               P1.reinit[xpred.rows,xpred.rows])
}
)

test_that("Reinitialised model with xpred, seasonal and AR1 uses prior information correctly", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  est.start.1 <- as.Date("2021-02-01")
  est.end.2 <- as.Date("2021-06-25")
  reinit.date <- as.Date("2021-04-21")
  q.default <- NULL
  sea.default <- 7
  
  data(gauteng_weather_2021, package = "tsgc")
  gauteng_weather <- gauteng_weather_2021[, c(1, 3)]
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default,
    xpred = gauteng_weather, start.date = est.start.1, end.date = reinit.date
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = cumulative_cases, q = q.default, sea.period = sea.default,
    xpred = gauteng_weather, start.date = est.start.1, end.date = est.end.2,
    reinit.date = reinit.date
  )
  
  i.reinit <- which(index(model_rei_base$Y)==reinit.date)
  lYy <- as.numeric(
    log(cumulative_cases[reinit.date]/diff(cumulative_cases)[reinit.date]))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit,]
  Tt <- est_base$output$model$T[,,1]
  Ptt <- est_base$output$Ptt[,,i.reinit]
  Rt <- est_base$output$model$R[,,1]
  Qt <- est_base$output$model$Q[,,1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base)=="level")
  slope.idx <- which(rownames(a1.base)=="slope")
  sea.rows <- grep("sea_trig",rownames(a1.base))
  xpred.rows <- grep("xpred",rownames(a1.base))
  choose <- c(level.idx, slope.idx, sea.rows, xpred.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[c(slope.idx,xpred.rows)] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.base[xpred.rows,xpred.rows] <- 0
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx,level.idx:slope.idx], 
               P1.reinit[level.idx:slope.idx,level.idx:slope.idx])
  expect_equal(P1.base[sea.rows,sea.rows],P1.reinit[sea.rows,sea.rows])
  expect_equal(P1.base[xpred.rows,xpred.rows],P1.reinit[xpred.rows,xpred.rows])
}
)

test_that("Model works with quarterly data", {
  data(nintendo_sales, package = "tsgc")
  wii <- nintendo_sales[, 1]
  
  est.start.q <- zoo::as.yearqtr("2006 Q4")
  est.end.q   <- zoo::as.yearqtr("2010 Q3")
  
  mod_wii <- tsgc::SSModelDynamicGompertz$new(
    Y = wii, sea.period = 4, start.date = est.start.q, end.date = est.end.q
  )
  res_wii <- mod_wii$estimate()
  
  expect_equal(length(res_wii$output$alphahat[16,]),
               2 + (res_wii$sea.period-1))
})

test_that("Model works with monthly data", {
  data(etrading_apps, package = "tsgc")
  Plus500 <- etrading_apps[, 1]

  est.start.m <- zoo::as.yearmon(2016)
  est.end.m   <- zoo::as.yearmon(2021)
  
  mod_500 <- tsgc::SSModelDynamicGompertz$new(
    Y = Plus500, sea.period = 12, start.date = est.start.m, end.date = est.end.m
  )
  res_500 <- mod_500$estimate()
  
  expect_equal(length(res_500$output$alphahat[16,]),
               2 + (res_500$sea.period-1))
})

test_that("Model works with annual data", {
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
  
  expect_equal(length(res_3ds$output$alphahat[6,]),2)
})
