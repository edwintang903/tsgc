library(KFAS)

test_that("tsgc gives same output as KFAS", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  dat <- get_timeframe(conv$series, conv$series$start, conv$series$start + 99)
  dat.ldl <- df2ldl(dat)
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = NULL,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  final.states.tsgc <- tsgc.est$output$alphahat
  
  kfas.model <- SSModel(as.matrix(idx_values(dat.ldl)) ~
                          SSMtrend(degree = 2,
                                   Q = list(matrix(0), matrix(NA))),
                        H = matrix(NA))
  kfas.fit <- fitSSM(kfas.model, inits = c(0, 0))
  kfas.out <- KFS(kfas.fit$model)
  final.states.kfas <- kfas.out$alphahat
  
  expect_equal(final.states.tsgc, final.states.kfas, tolerance = 1e-4)
})

test_that("tsgc enforces signal-to-noise restrictions correctly", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  dat <- get_timeframe(conv$series, conv$series$start, conv$series$start + 99)
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  
  q.out <- as.vector(tsgc.est$output$model$Q[2, 2, 1] / tsgc.est$output$model$H)
  
  expect_equal(q.out, q.choose)
})

test_that("tsgc gives same output as KFAS for fixed q", {
  updatesn <- function(pars, model, snr, order, index) {
    if (any(is.na(model$Q))) {
      Q <- as.matrix(model$Q[, , 1])
      naQd  <- which(is.na(diag(Q)))
      naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
      Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
      diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
      Q[naQd, naQd][naQnd] <- pars[length(naQd) + 1:length(naQnd)]
      model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
    }
    if (!identical(model$H, "Omitted") && any(is.na(model$H))) {
      H <- as.matrix(model$H[, , 1])
      naHd  <- which(is.na(diag(H)))
      naHnd <- which(upper.tri(H[naHd, naHd]) & is.na(H[naHd, naHd]))
      H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
      diag(H)[naHd] <-
        exp(0.5 * pars[length(naQd) + length(naQnd) + 1:length(naHd)])
      H[naHd, naHd][naHnd] <-
        pars[length(naQd) + length(naQnd) + length(naHd) + 1:length(naHnd)]
      model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
      model$Q[order, order, 1] <- snr * crossprod(H[index, index])
    }
    model
  }
  
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  dat <- get_timeframe(conv$series, conv$series$start, conv$series$start + 99)
  dat.ldl <- df2ldl(dat)
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = 0.005,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  final.states.tsgc <- tsgc.est$output$alphahat
  
  kfas.model <- SSModel(as.matrix(idx_values(dat.ldl)) ~
                          SSMtrend(degree = 2,
                                   Q = list(matrix(0), matrix(NA))),
                        H = matrix(NA))
  npar <- sum(is.na(kfas.model$Q)) + sum(is.na(kfas.model$H))
  update <- purrr::partial(updatesn, snr = 0.005, order = 2, index = 1)
  kfas.fit <- fitSSM(kfas.model, inits = rep(0, npar), updatefn = update)
  kfas.out <- KFS(kfas.fit$model)
  final.states.kfas <- kfas.out$alphahat
  
  expect_equal(final.states.tsgc, final.states.kfas, tolerance = 1e-4)
})

test_that("tsgc summary method functions", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  dat <- get_timeframe(conv$series, conv$series$start, conv$series$start + 99)
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  expect_no_error(expect_no_warning(tsgc.model$summary()))
})

test_that("tsgc print method functions", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  dat <- get_timeframe(conv$series, conv$series$start, conv$series$start + 99)
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  expect_no_error(expect_no_warning(tsgc.model$print()))
})


test_that("Model with seasonal has expected number of seasonal components", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  dat <- get_timeframe(conv$series, conv$series$start, conv$series$start + 99)
  
  sea.choose <- 7
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = 0.005,
                                           sea.period = sea.choose)
  tsgc.est <- tsgc.model$estimate()
  sea.est <- length(tsgc.est$output$model["a1", "seasonal"])
  
  expect_equal(sea.est, (sea.choose - 1))
})

test_that("Model with xpred has expected number of slope coefficients", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv_g <- xts_to_idx(gauteng[, 1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = idx_to_pos(conv_g$calendar, zoo::index(gauteng_weather_2021)[1]))
  
  est.start.1 <- idx_to_pos(conv_g$calendar, as.Date("2021-02-01"))
  est.end.1   <- idx_to_pos(conv_g$calendar, as.Date("2021-04-19"))
  
  model_weather <- tsgc::SSModelDynamicGompertz$new(
    Y = conv_g$series,
    xpred = conv_w$series,
    start = est.start.1,
    end   = est.end.1,
    sea.period = 0
  )
  
  est_weather <- model_weather$estimate()
  
  expect_equal(length(est_weather$output$model["a1", "regression"]),
               idx_ncol(conv_w$series))
})

test_that("Model with xpred and seasonal has expected number of
  seasonal and slope coefficients", {
    data(gauteng, package = "tsgc")
    data(gauteng_weather_2021, package = "tsgc")
    conv_g <- xts_to_idx(gauteng[, 1])
    conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = idx_to_pos(conv_g$calendar, zoo::index(gauteng_weather_2021)[1]))
    
    est.start.1 <- idx_to_pos(conv_g$calendar, as.Date("2021-02-01"))
    est.end.1   <- idx_to_pos(conv_g$calendar, as.Date("2021-04-19"))
    
    sea.choose <- 7
    
    model_weather <- tsgc::SSModelDynamicGompertz$new(
      Y = conv_g$series,
      xpred = conv_w$series,
      start = est.start.1,
      end   = est.end.1,
      sea.period = sea.choose
    )
    
    est_weather <- model_weather$estimate()
    
    expect_equal(length(est_weather$output$model["a1", "regression"]),
                 idx_ncol(conv_w$series))
    expect_equal(length(est_weather$output$model["a1", "seasonal"]),
                 (sea.choose - 1))
  })

test_that("Reinitialised model uses prior information correctly", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 0
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    start = est.start.1, end = reinit.pos
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    start = est.start.1, end = est.end.2,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  level.idx <- which(names(att) == "level")
  slope.idx <- which(names(att) == "slope")
  
  adj <- rep(0, length(att))
  adj[level.idx] <- lYy
  
  a1.base <- Tt %*% att + adj
  a1.base[slope.idx] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base, a1.reinit)
  expect_equal(P1.base, P1.reinit)
})

test_that("Reinitialised model with seasonal uses prior information correctly", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 7
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    start = est.start.1, end = reinit.pos
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    start = est.start.1, end = est.end.2,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base) == "level")
  slope.idx <- which(rownames(a1.base) == "slope")
  sea.rows <- grep("sea_trig", rownames(a1.base))
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[slope.idx] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base, a1.reinit)
  expect_equal(P1.base[level.idx:slope.idx, level.idx:slope.idx],
               P1.reinit[level.idx:slope.idx, level.idx:slope.idx])
  expect_equal(P1.base[sea.rows, sea.rows], P1.reinit[sea.rows, sea.rows])
})

test_that("Reinitialised model with xpred uses prior information correctly", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1]))
  
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 0
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    xpred = conv_w$series, start = est.start.1, end = reinit.pos
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    xpred = conv_w$series, start = est.start.1, end = est.end.2,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base) == "level")
  slope.idx <- which(rownames(a1.base) == "slope")
  xpred.rows <- grep("xpred", rownames(a1.base))
  choose <- c(level.idx, slope.idx, xpred.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[c(slope.idx, xpred.rows)] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.base[xpred.rows, xpred.rows] <- 0
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx, level.idx:slope.idx],
               P1.reinit[level.idx:slope.idx, level.idx:slope.idx])
  expect_equal(P1.base[xpred.rows, xpred.rows],
               P1.reinit[xpred.rows, xpred.rows])
})

test_that("Reinitialised model with xpred and seasonal uses prior information correctly", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1]))
  
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 7
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    xpred = conv_w$series, start = est.start.1, end = reinit.pos
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    xpred = conv_w$series, start = est.start.1, end = est.end.2,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base) == "level")
  slope.idx <- which(rownames(a1.base) == "slope")
  sea.rows <- grep("sea_trig", rownames(a1.base))
  xpred.rows <- grep("xpred", rownames(a1.base))
  choose <- c(level.idx, slope.idx, sea.rows, xpred.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[c(slope.idx, xpred.rows)] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.base[xpred.rows, xpred.rows] <- 0
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx, level.idx:slope.idx],
               P1.reinit[level.idx:slope.idx, level.idx:slope.idx])
  expect_equal(P1.base[sea.rows, sea.rows], P1.reinit[sea.rows, sea.rows])
  expect_equal(P1.base[xpred.rows, xpred.rows], P1.reinit[xpred.rows, xpred.rows])
})

test_that("Reinitialised model with AR1 (no xpred, no seasonal) uses prior information correctly", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 7
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    start = est.start.1, end = reinit.pos, ar1 = TRUE
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    start = est.start.1, end = est.end.2, ar1 = TRUE,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base) == "level")
  slope.idx <- which(rownames(a1.base) == "slope")
  sea.rows <- grep("sea_trig", rownames(a1.base))
  choose <- c(level.idx, slope.idx, sea.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[slope.idx] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx, level.idx:slope.idx],
               P1.reinit[level.idx:slope.idx, level.idx:slope.idx])
  expect_equal(P1.base[sea.rows, sea.rows], P1.reinit[sea.rows, sea.rows])
})

test_that("Reinitialised model with xpred and AR1 uses prior information correctly", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1]))
  
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 0
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default, ar1 = TRUE,
    xpred = conv_w$series, start = est.start.1, end = reinit.pos
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default, ar1 = TRUE,
    xpred = conv_w$series, start = est.start.1, end = est.end.2,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base) == "level")
  slope.idx <- which(rownames(a1.base) == "slope")
  xpred.rows <- grep("xpred", rownames(a1.base))
  choose <- c(level.idx, slope.idx, xpred.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[c(slope.idx, xpred.rows)] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.base[xpred.rows, xpred.rows] <- 0
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx, level.idx:slope.idx],
               P1.reinit[level.idx:slope.idx, level.idx:slope.idx])
  expect_equal(P1.base[xpred.rows, xpred.rows],
               P1.reinit[xpred.rows, xpred.rows])
})

test_that("Reinitialised model with xpred, seasonal and AR1 uses prior information correctly", {
  data(gauteng, package = "tsgc")
  data(gauteng_weather_2021, package = "tsgc")
  conv <- xts_to_idx(gauteng[, 1])
  conv_w <- xts_to_idx(gauteng_weather_2021[, c(1, 3)], start.pos = idx_to_pos(conv$calendar, zoo::index(gauteng_weather_2021)[1]))
  
  est.start.1 <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end.2   <- idx_to_pos(conv$calendar, as.Date("2021-06-25"))
  reinit.pos  <- idx_to_pos(conv$calendar, as.Date("2021-04-21"))
  q.default <- NULL
  sea.default <- 7
  
  model_rei_base <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    xpred = conv_w$series, start = est.start.1, end = reinit.pos
  )
  
  model_reinit <- tsgc::SSModelDynamicGompertz$new(
    Y = conv$series, q = q.default, sea.period = sea.default,
    xpred = conv_w$series, start = est.start.1, end = est.end.2,
    reinit.idx = reinit.pos
  )
  
  i.reinit <- which(idx_positions(model_rei_base$Y) == reinit.pos)
  lYy <- as.numeric(
    log(idx_values(conv$series[reinit.pos]) / idx_values(idx_diff(conv$series, 1L)[reinit.pos])))
  
  est_base <- model_rei_base$estimate()
  est_reinit <- model_reinit$estimate()
  
  att <- est_base$output$att[i.reinit, ]
  Tt <- est_base$output$model$T[, , 1]
  Ptt <- est_base$output$Ptt[, , i.reinit]
  Rt <- est_base$output$model$R[, , 1]
  Qt <- est_base$output$model$Q[, , 1]
  
  P1.base <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
  
  a1.base <- Tt %*% att
  level.idx <- which(rownames(a1.base) == "level")
  slope.idx <- which(rownames(a1.base) == "slope")
  sea.rows <- grep("sea_trig", rownames(a1.base))
  xpred.rows <- grep("xpred", rownames(a1.base))
  choose <- c(level.idx, slope.idx, sea.rows, xpred.rows)
  adj <- rep(0, length(a1.base))
  adj[level.idx] <- lYy
  a1.base <- a1.base + adj
  a1.base[c(slope.idx, xpred.rows)] <- 0
  a1.reinit <- est_reinit$output$model$a1
  
  P1.base[xpred.rows, xpred.rows] <- 0
  P1.reinit <- est_reinit$output$model$P1
  
  expect_equal(a1.base[choose], a1.reinit[choose])
  expect_equal(P1.base[level.idx:slope.idx, level.idx:slope.idx],
               P1.reinit[level.idx:slope.idx, level.idx:slope.idx])
  expect_equal(P1.base[sea.rows, sea.rows], P1.reinit[sea.rows, sea.rows])
  expect_equal(P1.base[xpred.rows, xpred.rows], P1.reinit[xpred.rows, xpred.rows])
})

test_that("Model works with quarterly, position-based data (no calendar attached)", {
  data(nintendo_sales, package = "tsgc")
  # yearqtr-indexed; build the calendar directly with unit = "quarters"
  # rather than via xts_to_idx() (which assumes a daily index).
  first_qtr <- zoo::index(nintendo_sales)[1]
  conv <- list(
    series = idx_series(zoo::coredata(nintendo_sales[, 1]), start = 1L),
    calendar = idx_calendar(anchor = zoo::as.Date(first_qtr), anchor_pos = 1L,
                            amount = 1, unit = "quarters", posixct = TRUE)
  )
  wii <- conv$series
  
  est.start.q <- idx_to_pos(conv$calendar, as.Date("2006-10-01"))
  est.end.q   <- idx_to_pos(conv$calendar, as.Date("2010-07-01"))
  
  mod_wii <- tsgc::SSModelDynamicGompertz$new(
    Y = wii, sea.period = 4, start = est.start.q, end = est.end.q
  )
  res_wii <- mod_wii$estimate()
  
  expect_equal(length(res_wii$output$alphahat[16, ]),
               2 + (res_wii$sea.period - 1))
})

test_that("Model works with monthly data converted via xts_to_idx", {
  data(etrading_apps, package = "tsgc")
  # yearmon-indexed; build the calendar directly with unit = "months".
  first_mon <- zoo::index(etrading_apps)[1]
  conv <- list(
    series = idx_series(zoo::coredata(etrading_apps[, 1]), start = 1L),
    calendar = idx_calendar(anchor = zoo::as.Date(first_mon), anchor_pos = 1L,
                            amount = 1, unit = "months", posixct = TRUE)
  )
  Plus500 <- conv$series
  
  est.start.m <- idx_to_pos(conv$calendar, zoo::as.Date(zoo::as.yearmon(2016)))
  est.end.m   <- idx_to_pos(conv$calendar, zoo::as.Date(zoo::as.yearmon(2021)))
  
  mod_500 <- tsgc::SSModelDynamicGompertz$new(
    Y = Plus500, sea.period = 12, start = est.start.m, end = est.end.m
  )
  res_500 <- mod_500$estimate()
  
  expect_equal(length(res_500$output$alphahat[16, ]),
               2 + (res_500$sea.period - 1))
})

test_that("Model works with a plain, non-calendar idx_series (arbitrary integer positions)", {
  set.seed(2)
  Y <- idx_series(cumsum(rpois(60, 6)) + 1, start = 1000L)
  
  mod <- tsgc::SSModelDynamicGompertz$new(Y = Y, sea.period = 0, q = 0.005,
                                          start = 1000L, end = 1050L)
  res <- mod$estimate()
  
  expect_true(inherits(res, "FilterResults"))
  expect_equal(res$index[1], 1000L)
  expect_equal(tail(res$index, 1), 1050L)
})

## Non-posixct calendar coverage: estimation equivalence, xpred.new ----

test_that("SSModelDynamicGompertz estimates identically regardless of calendar posixct flag", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  
  model_posix <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
                                            end = end.pos, calendar = conv$calendar)
  model_nonposix <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
                                               end = end.pos, calendar = cal_np)
  
  res_posix <- estimate(model_posix)
  res_nonposix <- estimate(model_nonposix)
  
  expect_false(res_nonposix$calendar$posixct)
  expect_equal(att(output(res_posix)), att(output(res_nonposix)))
  expect_equal(gety.hat(res_posix$predict_all(n.ahead = 7)),
               gety.hat(res_nonposix$predict_all(n.ahead = 7)))
})

test_that("SSModelDynamicGompertz works with a fully non-calendar (numeric-anchor, non-posixct) calendar", {
  set.seed(42)
  Y <- idx_series(cumsum(rpois(80, 8)) + 1, start = 1L)
  cal_ps <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                         unit = "picoseconds", posixct = FALSE)
  
  model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 60, calendar = cal_ps)
  res <- estimate(model)
  expect_s4_class(res, "FilterResults")
  expect_false(res$calendar$posixct)
  expect_no_error(res$predict_all(n.ahead = 7))
})

test_that("supplying res$xpred.new works when the fitted model used a non-posixct calendar", {
  set.seed(99)
  n <- 80
  y <- idx_series(cumsum(rpois(n, 8)) + 1, start = 1L)
  x <- idx_series(matrix(rnorm(n), ncol = 1), start = 1L)
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  
  model <- SSModelDynamicGompertz$new(Y = y, xpred = x, q = 0.005, end = 60, calendar = cal_np)
  res <- estimate(model)
  
  future_x <- idx_series(matrix(rnorm(10), ncol = 1), start = 61L)
  res$xpred.new <- future_x
  expect_no_error(res$predict_all(n.ahead = 10))
  expect_no_error(plot_forecast(res, n.ahead = 10))
})

test_that("estimation and predict_all give identical results whether calendar is a posixct=FALSE object or NULL entirely", {
  set.seed(25)
  Y <- idx_series(cumsum(rpois(60, 8)) + 1, start = 1L)
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  
  res_cal <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 50, calendar = cal_np)$estimate()
  res_nocal <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 50)$estimate()
  
  expect_equal(att(output(res_cal)), att(output(res_nocal)))
  expect_equal(res_cal$index, res_nocal$index)
  expect_equal(gety.hat(res_cal$predict_all(n.ahead = 7)),
               gety.hat(res_nocal$predict_all(n.ahead = 7)))
})