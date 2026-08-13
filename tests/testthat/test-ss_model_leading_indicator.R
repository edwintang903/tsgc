library(KFAS)

test_that("tsgc produces same LI output as KFAS", {
  set.seed(123)
  data(ukitaly, package = "tsgc")
  conv <- xts_to_idx(ukitaly)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2020-02-25"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2020-04-01"))
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = conv$series, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start = est.start,
                                      end = est.end, LeadIndCol = 1)
  tsgc_est <- tsgc_mod$estimate()
  
  y <- df2ldl_lead(conv$series, LeadIndCol = 1)
  
  y$newLead <- idx_lag(y$newLead, n.lag)
  y$LDLlead <- idx_lag(y$LDLlead, n.lag)
  y$cLead   <- idx_lag(y$cLead, n.lag)
  
  common_pos <- Reduce(intersect, lapply(y, idx_positions))
  combined_mat <- do.call(cbind, lapply(y, function(s) idx_values(s[common_pos])))
  colnames(combined_mat) <- names(y)
  y_combined <- idx_series(combined_mat, start = common_pos[1])
  
  finite_rows <- apply(idx_values(y_combined), 1, function(row) all(is.finite(row)))
  keep_pos <- idx_positions(y_combined)[finite_rows]
  y_clean <- y_combined[keep_pos]
  
  data_ldl <- get_timeframe(y_clean, est.start, est.end)
  data_mat <- idx_values(data_ldl)[, c("LDLlead", "LDLtarg")]
  
  kfas_mod <- SSModel(data_mat ~ SSMtrend(degree = 2,
                                          Q = matrix(c(0, 0, 0, NA), 2, 2),
                                          type = "common") +
                        SSMtrend(degree = 1, Q = matrix(NA), index = 1),
                      H = matrix(c(NA, 0, 0, NA), 2, 2))
  npar <- sum(is.na(kfas_mod$Q)) + sum(is.na(kfas_mod$H))
  kfas_fit <- fitSSM(kfas_mod, rep(0, npar))
  kfas_est <- KFS(kfas_fit$model)
  
  expect_equal(unname(as.matrix(tsgc_est$output$alphahat)),
               unname(as.matrix(kfas_est$alphahat)))
})

test_that("tsgc produces same LI output as KFAS with fixed q", {
  data(ukitaly, package = "tsgc")
  conv <- xts_to_idx(ukitaly)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2020-02-25"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2020-04-01"))
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = conv$series, n.lag = n.lag, q = 0.005,
                                      sea.period = 0, start = est.start,
                                      end = est.end, LeadIndCol = 1)
  tsgc_est <- tsgc_mod$estimate()
  
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
  updateli <- updatesn %>% purrr::partial(snr = 0.005, order = 2, index = 2)
  
  y <- df2ldl_lead(conv$series, LeadIndCol = 1)
  y$newLead <- idx_lag(y$newLead, n.lag)
  y$LDLlead <- idx_lag(y$LDLlead, n.lag)
  y$cLead   <- idx_lag(y$cLead, n.lag)
  
  common_pos <- Reduce(intersect, lapply(y, idx_positions))
  combined_mat <- do.call(cbind, lapply(y, function(s) idx_values(s[common_pos])))
  colnames(combined_mat) <- names(y)
  y_combined <- idx_series(combined_mat, start = common_pos[1])
  
  finite_rows <- apply(idx_values(y_combined), 1, function(row) all(is.finite(row)))
  keep_pos <- idx_positions(y_combined)[finite_rows]
  y_clean <- y_combined[keep_pos]
  
  data_ldl <- get_timeframe(y_clean, est.start, est.end)
  data_mat <- idx_values(data_ldl)[, c("LDLlead", "LDLtarg")]
  
  kfas_mod <- SSModel(data_mat ~ SSMtrend(degree = 2,
                                          Q = matrix(c(0, 0, 0, NA), 2, 2),
                                          type = "common") +
                        SSMtrend(degree = 1, Q = matrix(NA), index = 1),
                      H = matrix(c(NA, 0, 0, NA), 2, 2))
  npar <- sum(is.na(kfas_mod$Q)) + sum(is.na(kfas_mod$H))
  kfas_fit <- fitSSM(kfas_mod, rep(0, npar), updatefn = updateli)
  kfas_est <- KFS(kfas_fit$model)
  
  tsgc_snr <- tsgc_est$output$model$Q[2, 2, 1] / tsgc_est$output$model$H[2, 2, 1]
  kfas_snr <- kfas_est$model$Q[2, 2, 1] / kfas_est$model$H[2, 2, 1]
  
  expect_equal(tsgc_snr, 0.005)
  expect_equal(kfas_snr, 0.005)
  expect_equal(unname(as.matrix(tsgc_est$output$alphahat)),
               unname(as.matrix(kfas_est$alphahat)))
})

test_that("Summary method works", {
  data(ukitaly, package = "tsgc")
  conv <- xts_to_idx(ukitaly)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2020-02-25"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2020-04-01"))
  
  tsgc_mod <- SSModelLeadingIndicator(Y = conv$series, n.lag = 14, q = NULL,
                                      sea.period = 0, start = est.start,
                                      end = est.end, LeadIndCol = 1)
  
  expect_no_error(expect_no_warning(tsgc_mod$summary()))
})

test_that("Print method works", {
  data(ukitaly, package = "tsgc")
  conv <- xts_to_idx(ukitaly)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2020-02-25"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2020-04-01"))
  
  tsgc_mod <- SSModelLeadingIndicator(Y = conv$series, n.lag = 14, q = NULL,
                                      sea.period = 0, start = est.start,
                                      end = est.end, LeadIndCol = 1)
  
  expect_no_error(expect_no_warning(tsgc_mod$print()))
})

test_that("Leading indicator model + seasonal has correct number of elements", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  
  est.start.eng <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  est.end.eng   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  sea <- 7
  
  mod <- SSModelLeadingIndicator$new(conv$series, n.lag = 5, start = est.start.eng,
                                     end = est.end.eng, sea.period = sea)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat), 3 + 2 * (sea - 1))
})

test_that("LI model + xpred_lead + seasonal has correct number of elements", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  conv_xp <- xts_to_idx(england_weather_2021[, 1:4], start.pos = idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1]))
  
  est.start.eng <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  est.end.eng   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  sea <- 7
  
  mod <- SSModelLeadingIndicator$new(conv$series, n.lag = 5, start = est.start.eng,
                                     end = est.end.eng, sea.period = sea,
                                     xpred_lead = conv_xp$series)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat), 3 + 2 * (sea - 1) + idx_ncol(conv_xp$series))
})

test_that("LI model + xpred_targ has correct number of elements", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  conv_xp <- xts_to_idx(england_weather_2021[, 1:4], start.pos = idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1]))
  
  est.start.eng <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  est.end.eng   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  sea <- 0
  
  mod <- SSModelLeadingIndicator$new(conv$series, n.lag = 5, start = est.start.eng,
                                     end = est.end.eng, sea.period = sea,
                                     xpred_targ = conv_xp$series)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat), 3 + idx_ncol(conv_xp$series))
})

test_that("LI + xpred_lead + xpred_targ + seasonal has correct number of elements", {
  data(england, package = "tsgc")
  data(england_weather_2021, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  conv_xp <- xts_to_idx(england_weather_2021[, 1:4], start.pos = idx_to_pos(conv$calendar, zoo::index(england_weather_2021)[1]))
  
  est.start.eng <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  est.end.eng   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  sea <- 7
  
  mod <- SSModelLeadingIndicator$new(conv$series, n.lag = 5, start = est.start.eng,
                                     end = est.end.eng, sea.period = sea,
                                     xpred_lead = conv_xp$series, xpred_targ = conv_xp$series)
  res <- mod$estimate()
  expect_equal(ncol(res$output$alphahat), 3 + 2 * (sea - 1) + 2 * idx_ncol(conv_xp$series))
})

test_that("LI with quarterly data has correct number of components", {
  data(nintendo_sales, package = "tsgc")
  first_qtr <- zoo::index(nintendo_sales)[1]
  conv <- list(
    series = idx_series(zoo::coredata(nintendo_sales[, c("wii", "switch_all")]), start = 1L),
    calendar = idx_calendar(anchor = zoo::as.Date(first_qtr), anchor_pos = 1L,
                            amount = 1, unit = "quarters", posixct = TRUE)
  )
  
  sea <- 4
  est.start.q2 <- idx_to_pos(conv$calendar, as.Date("2017-01-01"))
  est.end.q2   <- idx_to_pos(conv$calendar, as.Date("2019-10-01"))
  n.lag.q      <- est.start.q2 - idx_to_pos(conv$calendar, as.Date("2006-10-01"))
  
  mod_switch <- tsgc::SSModelLeadingIndicator$new(
    Y = conv$series, sea.period = sea, n.lag = n.lag.q,
    start = est.start.q2, end = est.end.q2
  )
  res <- mod_switch$estimate()
  
  expect_equal(ncol(res$output$alphahat), 3 + 2 * (sea - 1))
})

test_that("LI works with a plain, non-calendar idx_series (arbitrary integer positions)", {
  set.seed(4)
  lead <- cumsum(rpois(150, 6)) + 1
  targ <- cumsum(rpois(150, 8)) + 1
  Y <- idx_series(cbind(lead, targ), start = 1L)
  
  mod <- SSModelLeadingIndicator(Y = Y, n.lag = 5, q = NULL, LeadIndCol = 1,
                                 sea.period = 0, start = 1L, end = 100L)
  res <- estimate(mod)
  
  expect_true(inherits(res, "FilterResultsLI"))
  expect_true(res$start > 1L)
  expect_equal(res$end, 100L)
})

## Non-posixct calendar coverage: estimation equivalence, quarterly/monthly frequency ----

test_that("SSModelLeadingIndicator estimates identically regardless of calendar posixct flag", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england[, 1:2])
  
  cal_np <- idx_calendar(anchor = zoo::index(england)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  
  start.pos <- idx_to_pos(conv$calendar, as.Date("2021-04-30"))
  end.pos   <- idx_to_pos(conv$calendar, as.Date("2021-07-24"))
  
  model_posix <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 7,
    start = start.pos, end = end.pos, calendar = conv$calendar
  )
  model_nonposix <- SSModelLeadingIndicator$new(
    conv$series, n.lag = 4, sea.period = 7,
    start = start.pos, end = end.pos, calendar = cal_np
  )
  
  res_posix <- estimate(model_posix)
  res_nonposix <- estimate(model_nonposix)
  
  expect_false(res_nonposix$calendar$posixct)
  expect_equal(att(output(res_posix)), att(output(res_nonposix)))
})

test_that("quarterly-frequency data works under a non-posixct calendar (plain integer step, not calendar-aware quarter stepping)", {
  data(nintendo_sales, package = "tsgc")
  y_q_xts <- nintendo_sales[, c("wii", "switch_all")]
  nintendo_idx <- idx_series(zoo::coredata(y_q_xts), start = 1L)
  
  nintendo_cal_np <- idx_calendar(
    anchor = zoo::as.Date(zoo::index(y_q_xts)[1]),
    anchor_pos = 1L, amount = 1, unit = "quarters", posixct = FALSE
  )
  expect_error(idx_to_pos(nintendo_cal_np, "2017-01-01"), "not a calendar-anchored")
  
  nintendo_cal_helper <- idx_calendar(
    anchor = zoo::as.Date(zoo::index(y_q_xts)[1]),
    anchor_pos = 1L, amount = 1, unit = "quarters", posixct = TRUE
  )
  n.lag.q <- idx_to_pos(nintendo_cal_helper, "2017-01-01") - idx_to_pos(nintendo_cal_helper, "2006-10-01")
  start.q <- idx_to_pos(nintendo_cal_helper, "2017-01-01")
  end.q   <- idx_to_pos(nintendo_cal_helper, "2019-10-01")
  
  mod.q <- SSModelLeadingIndicator$new(
    Y = nintendo_idx, sea.period = 4, n.lag = n.lag.q,
    start = start.q, end = end.q, calendar = nintendo_cal_np
  )
  res.q <- estimate(mod.q)
  expect_no_error(tsgc::plot_log_forecast(res.q, Y = nintendo_idx, n.ahead = 8))
  expect_no_error(plot_gy_ci(res.q, axis = idx_axis_opts(mode = "steps")))
})

test_that("monthly-frequency data works under a non-posixct calendar", {
  data(etrading_apps, package = "tsgc")
  y_m_xts <- etrading_apps[, c("DEGIRO", "AvaTrade")]
  etrading_idx <- idx_series(zoo::coredata(y_m_xts), start = 1L)
  
  etrading_cal_np <- idx_calendar(
    anchor = zoo::as.Date(zoo::index(y_m_xts)[1]),
    anchor_pos = 1L, amount = 1, unit = "months", posixct = FALSE
  )
  etrading_cal_helper <- idx_calendar(
    anchor = zoo::as.Date(zoo::index(y_m_xts)[1]),
    anchor_pos = 1L, amount = 1, unit = "months", posixct = TRUE
  )
  n.lag.m <- idx_to_pos(etrading_cal_helper, "2017-07-01") - idx_to_pos(etrading_cal_helper, "2017-01-01")
  start.m <- idx_to_pos(etrading_cal_helper, "2017-07-01")
  end.m   <- idx_to_pos(etrading_cal_helper, "2021-02-01")
  
  mod.m <- SSModelLeadingIndicator$new(
    Y = etrading_idx, sea.period = 12, n.lag = n.lag.m,
    start = start.m, end = end.m, calendar = etrading_cal_np
  )
  res.m <- estimate(mod.m)
  expect_no_error(tsgc::plot_forecast(res.m, n.ahead = 4))
  expect_no_error(plot_forecast(res.m, n.ahead = 4, axis = idx_axis_opts(mode = "time_since")))
})

test_that("repeated LI estimation does not shift or truncate lead regressors", {
  set.seed(2)
  Y <- idx_series(cbind(cumsum(rpois(150, 6)) + 1,
                        cumsum(rpois(150, 8)) + 1))
  xpred <- idx_series(matrix(rnorm(300), ncol = 2))
  mod <- SSModelLeadingIndicator(
    Y, n.lag = 5, sea.period = 0, start = 20, end = 100,
    xpred_lead = xpred
  )
  original_xpred <- mod$xpred_lead
  
  first <- estimate(mod)
  second <- estimate(mod)
  
  expect_identical(mod$xpred_lead, original_xpred)
  expect_equal(first$output$alphahat, second$output$alphahat)
})
