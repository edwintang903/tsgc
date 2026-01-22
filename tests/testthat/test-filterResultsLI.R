library(testthat)
library(tsgc)
library(KFAS)
library(xts)
library(zoo)

data(england, package = "tsgc")
Y_eng <- england[, 1:2]

est_start <- as.Date("2021-04-30")
est_end   <- as.Date("2021-07-24")
n_lag     <- 4

make_england_res <- function(sea.period, start.date, end.date, n.lag = n_lag) {
  mod <- SSModelLeadingIndicator(
    Y          = Y_eng,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = sea.period,
    start.date = start.date,
    end.date   = end.date
  )
  estimate(mod)
}

# Base objects reused across tests (avoid re-estimation in each test_that)
res_base  <- make_england_res(sea.period = 7, start.date = est_start, end.date = est_end)
res_nosea <- make_england_res(sea.period = 0, start.date = est_start, end.date = est_end)

# A model estimated near the end of the dataset: remaining_data < n.lag and short holdout
last_date <- tail(index(Y_eng), 1)
est_end_short   <- last_date - 2
est_start_short <- est_end_short - 60
res_short_holdout <- make_england_res(sea.period = 7, start.date = est_start_short, end.date = est_end_short)

# ------------------------------------------------------------------------------
# Synthetic builder for branch coverage:
# - xpred branches in predict_all() + plot_log_forecast()
# - non-daily resolution branches in predict_level(), summary(), plots
# ------------------------------------------------------------------------------
make_synthetic_FilterResultsLI <- function(
    resolution = c("daily", "monthly", "quarterly", "yearly"),
    sea.period = 0,
    xpred_logical = c(FALSE, FALSE),
    n.lag = 4,
    n_total = 80,
    n_est = 50
) {
  resolution <- match.arg(resolution)
  stopifnot(length(xpred_logical) == 2)
  
  idx <- switch(
    resolution,
    daily     = seq.Date(as.Date("2020-01-01"), by = "day", length.out = n_total),
    monthly   = as.yearmon(seq.Date(as.Date("2020-01-01"), by = "month", length.out = n_total)),
    quarterly = as.yearqtr(seq.Date(as.Date("2020-01-01"), by = "3 months", length.out = n_total)),
    yearly    = as.yearmon(seq.Date(as.Date("2000-01-01"), by = "year", length.out = n_total))
  )
  
  # Deterministic, strictly positive cumulated series
  lead <- 1000 + cumsum(10 + (seq_len(n_total) %% 3))
  targ <- 2000 + cumsum(6  + (seq_len(n_total) %% 2))
  Y <- xts(cbind(lead, targ), order.by = idx)
  colnames(Y) <- c("lead", "targ")
  
  # Construct the fields FilterResultsLI expects: cLead/cTarg/newTarg/LDLlead/LDLtarg
  newLead <- c(NA_real_, diff(lead))
  newTarg <- c(NA_real_, diff(targ))
  LDLlead <- c(NA_real_, diff(log(lead)))
  LDLtarg <- c(NA_real_, diff(log(targ)))
  
  data_xts <- xts(
    cbind(
      cLead   = as.numeric(lead),
      cTarg   = as.numeric(targ),
      newLead = as.numeric(newLead),
      newTarg = as.numeric(newTarg),
      LDLlead = as.numeric(LDLlead),
      LDLtarg = as.numeric(LDLtarg)
    ),
    order.by = idx
  )
  
  # Estimation window avoids the first NA row created by diff()
  stopifnot(n_est + 1 < n_total)
  start.date <- idx[2]
  end.date   <- idx[n_est + 1]
  y_est <- as.matrix(data_xts[2:(n_est + 1), c("LDLlead", "LDLtarg")])
  
  # Exogenous predictors (aligned to full sample)
  x_lead_full <- xts(matrix(seq_len(n_total), ncol = 1), order.by = idx)
  x_targ_full <- xts(matrix(100 + seq_len(n_total), ncol = 1), order.by = idx)
  colnames(x_lead_full) <- "x_lead"
  colnames(x_targ_full) <- "x_targ"
  x_lead_est <- as.matrix(x_lead_full[2:(n_est + 1)])
  x_targ_est <- as.matrix(x_targ_full[2:(n_est + 1)])
  colnames(x_lead_est) <- "x_lead"
  colnames(x_targ_est) <- "x_targ"
  
  # Fixed-parameter SSModel compatible with FilterResultsLI$predict_all()
  Q_slope  <- 0.01
  Q_trend1 <- 0.01
  H_mat    <- diag(c(0.1, 0.1))
  
  comps <- SSMtrend(
    degree = 2,
    Q      = matrix(c(0, 0, 0, Q_slope), 2, 2),
    type   = "common"
  )
  
  if (sea.period >= 2) {
    comps <- comps + SSMseasonal(
      sea.period,
      Q        = matrix(0, 2, 2),
      sea.type = "trigonometric",
      type     = "distinct"
    )
  }
  
  comps <- comps + SSMtrend(degree = 1, Q = matrix(Q_trend1), index = 1)
  
  if (isTRUE(xpred_logical[1])) {
    comps <- comps + SSMregression(~x_lead_est, type = "distinct", index = 1)
  }
  if (isTRUE(xpred_logical[2])) {
    comps <- comps + SSMregression(~x_targ_est, type = "distinct", index = 2)
  }
  
  model <- SSModel(y_est ~ comps, H = H_mat)
  out   <- KFS(model, smoothing = c("state", "signal"))
  
  res <- FilterResultsLI$new(
    data_xts       = data_xts,
    output         = out,
    n.lag          = n.lag,
    sea.period     = sea.period,
    LeadIndCol     = 1,
    xpred_logical  = as.logical(xpred_logical),
    start.date     = start.date,
    end.date       = end.date,
    xpred_lead.new = if (isTRUE(xpred_logical[1])) x_lead_full else NULL,
    xpred_targ.new = if (isTRUE(xpred_logical[2])) x_targ_full else NULL
  )
  
  # Force desired resolution label to exercise the resolution-specific branches
  res$resolution <- resolution
  
  invisible(list(res = res, Y = Y))
}

# ------------------------------------------------------------------------------
# predict_level: baseline test (your existing manual construction check)
# ------------------------------------------------------------------------------
test_that("FilterResultsLI::predict_level matches manual construction (daily, seasonal, no xpred)", {
  nf <- 7
  
  kfas_fc <- res_base$predict_all(
    n.ahead          = nf,
    sea.on           = TRUE,
    return.all       = FALSE,
    confidence.level = 0.68
  )$y.hat.kfas
  
  delta_fit <- as.vector(kfas_fc$LDLtarg[, 1])
  delta_lwr <- as.vector(kfas_fc$LDLtarg[, 2])
  delta_upr <- as.vector(kfas_fc$LDLtarg[, 3])
  
  last_row <- get_timeframe(res_base$data_xts, res_base$end.date)[1, ]
  YT       <- as.numeric(last_row$cTarg)
  
  cp_fit   <- cumprod(1 + exp(delta_fit[1:(nf - 1)]))
  mult_fit <- c(1, cp_fit)
  forc_fit <- YT * exp(delta_fit) * mult_fit
  
  cp_lwr   <- cumprod(1 + exp(delta_lwr[1:(nf - 1)]))
  mult_lwr <- c(1, cp_lwr)
  forc_lwr <- YT * exp(delta_lwr) * mult_lwr
  
  cp_upr   <- cumprod(1 + exp(delta_upr[1:(nf - 1)]))
  mult_upr <- c(1, cp_upr)
  forc_upr <- YT * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res_base$predict_level(n.ahead = nf, confidence.level = 0.68, sea.on = TRUE)
  
  expect_equal(as.numeric(forc_tsgc$forc), round(forc_fit, 2))
  expect_equal(as.numeric(forc_tsgc$lwr),  round(forc_lwr, 2))
  expect_equal(as.numeric(forc_tsgc$upr),  round(forc_upr, 2))
})

test_that("FilterResultsLI::predict_level handles n.ahead == 1 (unity branch)", {
  fc1 <- res_base$predict_level(n.ahead = 1, sea.on = TRUE)
  
  expect_true(xts::is.xts(fc1))
  expect_equal(NROW(fc1), 1)
  expect_equal(colnames(fc1), c("forc", "lwr", "upr"))
  
  expected_first <- seq_dates(res_base$end.date, length.out = 2, resolution = res_base$resolution)[2]
  expect_equal(index(fc1)[1], expected_first)
})

test_that("FilterResultsLI::predict_level matches manual construction (daily, trend-only / sea.on=FALSE)", {
  nf <- 7
  cl <- 0.50  # non-default CI level exercises confidence.level plumbing for sea.on=FALSE
  
  kfas_fc <- res_base$predict_all(
    n.ahead          = nf,
    sea.on           = FALSE,
    return.all       = FALSE,
    confidence.level = cl
  )$y.hat.kfas
  
  delta_fit <- as.vector(kfas_fc$LDLtarg[, 1])
  delta_lwr <- as.vector(kfas_fc$LDLtarg[, 2])
  delta_upr <- as.vector(kfas_fc$LDLtarg[, 3])
  
  last_row <- get_timeframe(res_base$data_xts, res_base$end.date)[1, ]
  YT       <- as.numeric(last_row$cTarg)
  
  cp_fit   <- cumprod(1 + exp(delta_fit[1:(nf - 1)]))
  mult_fit <- c(1, cp_fit)
  forc_fit <- YT * exp(delta_fit) * mult_fit
  
  cp_lwr   <- cumprod(1 + exp(delta_lwr[1:(nf - 1)]))
  mult_lwr <- c(1, cp_lwr)
  forc_lwr <- YT * exp(delta_lwr) * mult_lwr
  
  cp_upr   <- cumprod(1 + exp(delta_upr[1:(nf - 1)]))
  mult_upr <- c(1, cp_upr)
  forc_upr <- YT * exp(delta_upr) * mult_upr
  
  forc_tsgc <- res_base$predict_level(n.ahead = nf, confidence.level = cl, sea.on = FALSE)
  
  expect_equal(as.numeric(forc_tsgc$forc), round(forc_fit, 2))
  expect_equal(as.numeric(forc_tsgc$lwr),  round(forc_lwr, 2))
  expect_equal(as.numeric(forc_tsgc$upr),  round(forc_upr, 2))
})

# ------------------------------------------------------------------------------
# predict_all: return.all, sea.on, n.ahead edge, sea.period==0, remaining_data edge
# ------------------------------------------------------------------------------
test_that("FilterResultsLI::predict_all respects return.all and sea.on flags", {
  nf <- 7
  
  out_all <- res_base$predict_all(n.ahead = nf, sea.on = TRUE,  return.all = TRUE)
  out_fc  <- res_base$predict_all(n.ahead = nf, sea.on = TRUE,  return.all = FALSE)
  out_tr  <- res_base$predict_all(n.ahead = nf, sea.on = FALSE, return.all = FALSE)
  
  expect_true(xts::is.xts(out_all$y.hat))
  expect_true(xts::is.xts(out_fc$y.hat))
  expect_true(xts::is.xts(out_tr$y.hat))
  
  n_est <- attr(res_base$output$model, "n")
  expect_equal(NROW(out_all$y.hat), n_est + nf)
  expect_equal(NROW(out_fc$y.hat), nf)
  expect_equal(NROW(out_tr$y.hat), nf)
  
  expect_equal(index(out_all$y.hat)[1], res_base$start.date)
  
  expect_true(is.matrix(out_fc$a.t.t))
  expect_true(is.array(out_fc$P.t.t))
})

test_that("FilterResultsLI::predict_all handles n.ahead <= n.lag (future_rows = n.ahead)", {
  out2 <- res_base$predict_all(n.ahead = 2, sea.on = TRUE, return.all = FALSE)
  expect_equal(NROW(out2$y.hat), 2)
})

test_that("FilterResultsLI::predict_all works when sea.period == 0 (no seasonal component)", {
  nf <- 7
  out <- res_nosea$predict_all(n.ahead = nf, sea.on = TRUE, return.all = FALSE)
  expect_equal(NROW(out$y.hat), nf)
  expect_no_error(res_nosea$predict_level(n.ahead = nf, sea.on = TRUE))
})

test_that("FilterResultsLI::predict_all works when remaining_data < min(n.ahead, n.lag)", {
  nf <- 7
  out <- res_short_holdout$predict_all(n.ahead = nf, sea.on = TRUE, return.all = FALSE)
  expect_equal(NROW(out$y.hat), nf)
})

# ------------------------------------------------------------------------------
# get_growth_y + get_gy_ci: all (smoothed/filtered) cases + confidence.level
# ------------------------------------------------------------------------------
test_that("FilterResultsLI::get_growth_y toggles smoothed and return.components (all cases)", {
  smooth_all  <- res_base$get_growth_y(smoothed = TRUE,  return.components = TRUE)
  smooth_only <- res_base$get_growth_y(smoothed = TRUE,  return.components = FALSE)
  filt_only   <- res_base$get_growth_y(smoothed = FALSE, return.components = FALSE)
  filt_all    <- res_base$get_growth_y(smoothed = FALSE, return.components = TRUE)
  
  expect_equal(length(smooth_all), 3)
  expect_true(xts::is.xts(smooth_only))
  expect_false(is.list(filt_only))
  
  expect_equal(names(smooth_all[[1]]), "smoothed gy.t")
  expect_equal(names(smooth_all[[2]]), "smoothed g.t")
  expect_equal(names(smooth_all[[3]]), "smoothed gamma.t")
  
  expect_equal(names(filt_all[[1]]), "filtered gy.t")
  expect_equal(names(filt_all[[2]]), "filtered g.t")
  expect_equal(names(filt_all[[3]]), "filtered gamma.t")
})

test_that("FilterResultsLI::get_gy_ci respects confidence.level and smoothed flag", {
  ci_68_f <- res_base$get_gy_ci(smoothed = FALSE, confidence.level = 0.68)
  ci_95_f <- res_base$get_gy_ci(smoothed = FALSE, confidence.level = 0.95)
  ci_68_s <- res_base$get_gy_ci(smoothed = TRUE,  confidence.level = 0.68)
  
  expect_equal(colnames(ci_68_f), c("fit", "lower", "upper"))
  expect_false(isTRUE(all.equal(ci_68_f, ci_95_f)))
  expect_false(isTRUE(all.equal(ci_68_f, ci_68_s)))
})

# ------------------------------------------------------------------------------
# print/summary + plots: optional arguments + non-default branches
# ------------------------------------------------------------------------------
test_that("FilterResultsLI print/summary/print_estimation_results do not error", {
  expect_no_error(res_base$print_estimation_results())
  expect_no_error(print(res_base))
  expect_no_error(summary(res_base))
})

test_that("FilterResultsLI plotting methods cover optional arguments and smoothed variants", {
  expect_no_error(res_base$plot_forecast(plt.start.date = res_base$end.date - 14))
  expect_no_error(res_base$plot_log_forecast(Y = Y_eng, plt.start.date = res_base$end.date - 14, caption = ""))
  expect_no_error(res_base$plot_gy_components(plt.start.date = res_base$end.date - 14, smoothed = TRUE))
  expect_no_error(res_base$plot_gy_ci(plt.start.date = res_base$end.date - 14, smoothed = TRUE,
                                      series.name = "target", pad.right = 7))
})

# ------------------------------------------------------------------------------
# plot_holdout + mapes: error + warning branches (zeros), plus bounds on coverage
# ------------------------------------------------------------------------------
test_that("FilterResultsLI::plot_holdout errors when holdout shorter than n.ahead", {
  expect_error(
    res_short_holdout$plot_holdout(Y = Y_eng, n.ahead = 10),
    "shorter than n.ahead"
  )
})

test_that("FilterResultsLI::plot_holdout warns when validation data contains zeros", {
  Y_zero <- Y_eng
  end_pos <- which.max(index(Y_eng) == res_base$end.date)
  holdout_date <- index(Y_eng)[end_pos + 1]
  Y_zero[holdout_date, 2] <- Y_zero[res_base$end.date, 2]  # flatten cumulative => newTarg=0
  
  expect_warning(
    p <- res_base$plot_holdout(Y = Y_zero, n.ahead = 1),
    "contains zeros"
  )
  expect_true(inherits(p, "ggplot"))
})

test_that("FilterResultsLI::mapes returns 5 metrics; warns on zeros; coverage is bounded", {
  errs <- res_base$mapes(n.ahead = 7, Y = Y_eng)
  
  expect_true(all(c("mape", "smape", "mae", "rmse", "coverage") %in% names(errs)))
  expect_equal(length(errs), 5)
  expect_true(is.numeric(errs$coverage))
  expect_gte(errs$coverage, 0)
  expect_lte(errs$coverage, 100)
  
  Y_zero <- Y_eng
  end_pos <- which.max(index(Y_eng) == res_base$end.date)
  holdout_date <- index(Y_eng)[end_pos + 1]
  Y_zero[holdout_date, 2] <- Y_zero[res_base$end.date, 2]
  expect_warning(res_base$mapes(n.ahead = 1, Y = Y_zero), "contains zeros")
})

test_that("FilterResultsLI::plot_holdout works in a normal (non-error) case", {
  expect_no_error(res_base$plot_holdout(Y = Y_eng, n.ahead = 7))
})

# ------------------------------------------------------------------------------
# Non-daily resolutions: predict_level date increments + summary/plot date branches
# ------------------------------------------------------------------------------
test_that("Non-daily resolutions: predict_level increments + plot/summary branches", {
  # Quarterly
  syn_qtr <- make_synthetic_FilterResultsLI(resolution = "quarterly", sea.period = 0, xpred_logical = c(FALSE, FALSE),
                                            n_total = 40, n_est = 20)
  res_qtr <- syn_qtr$res
  fc_qtr  <- res_qtr$predict_level(n.ahead = 3, sea.on = TRUE)
  expected_qtr_first <- seq_dates(res_qtr$end.date, length.out = 2, resolution = res_qtr$resolution)[2]
  expect_equal(index(fc_qtr)[1], expected_qtr_first)
  expect_no_error(summary(res_qtr))
  expect_no_error(res_qtr$plot_forecast(n.ahead = 3))
  expect_no_error(res_qtr$plot_gy_components())
  expect_no_error(res_qtr$plot_gy_ci())
  
  # Monthly
  syn_mon <- make_synthetic_FilterResultsLI(resolution = "monthly", sea.period = 0, xpred_logical = c(FALSE, FALSE),
                                            n_total = 36, n_est = 18)
  res_mon <- syn_mon$res
  fc_mon  <- res_mon$predict_level(n.ahead = 3, sea.on = TRUE)
  expected_mon_first <- seq_dates(res_mon$end.date, length.out = 2, resolution = res_mon$resolution)[2]
  expect_equal(index(fc_mon)[1], expected_mon_first)
  expect_no_error(summary(res_mon))
  expect_no_error(res_mon$plot_forecast(n.ahead = 3))
  expect_no_error(res_mon$plot_gy_components())
  expect_no_error(res_mon$plot_gy_ci())
  
  # Yearly
  syn_yr <- make_synthetic_FilterResultsLI(resolution = "yearly", sea.period = 0, xpred_logical = c(FALSE, FALSE),
                                           n_total = 25, n_est = 12)
  res_yr <- syn_yr$res
  fc_yr  <- res_yr$predict_level(n.ahead = 3, sea.on = TRUE)
  expected_yr_first <- seq_dates(res_yr$end.date, length.out = 2, resolution = res_yr$resolution)[2]
  expect_equal(index(fc_yr)[1], expected_yr_first)
  expect_no_error(summary(res_yr))
  expect_no_error(res_yr$plot_forecast(n.ahead = 3))
  expect_no_error(res_yr$plot_gy_components())
  expect_no_error(res_yr$plot_gy_ci())
})

# ------------------------------------------------------------------------------
# Xpred: error branches + all (sea.period <2 / >=2) x (lead only/targ only/both)
# ------------------------------------------------------------------------------
test_that("predict_all xpred branches: errors when xpred_*.new missing", {
  syn_lead <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 0, xpred_logical = c(TRUE, FALSE),
                                             n_total = 80, n_est = 50)
  res_lead <- syn_lead$res
  res_lead$xpred_lead.new <- NULL
  expect_error(res_lead$predict_all(n.ahead = 3), "xpred_lead.new not provided")
  
  syn_targ <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 0, xpred_logical = c(FALSE, TRUE),
                                             n_total = 80, n_est = 50)
  res_targ <- syn_targ$res
  res_targ$xpred_targ.new <- NULL
  expect_error(res_targ$predict_all(n.ahead = 3), "xpred_targ.new not provided")
})

test_that("predict_all xpred branches: lead-only, targ-only, both; with and without seasonality", {
  nf <- 3
  
  # sea.period < 2 : no seasonality branch in xpred path
  syn0_lead <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 0, xpred_logical = c(TRUE, FALSE),
                                              n_total = 80, n_est = 50)
  syn0_targ <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 0, xpred_logical = c(FALSE, TRUE),
                                              n_total = 80, n_est = 50)
  syn0_both <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 0, xpred_logical = c(TRUE, TRUE),
                                              n_total = 80, n_est = 50)
  
  out0_lead <- syn0_lead$res$predict_all(n.ahead = nf, sea.on = TRUE,  return.all = FALSE)
  out0_targ <- syn0_targ$res$predict_all(n.ahead = nf, sea.on = TRUE,  return.all = FALSE)
  out0_both <- syn0_both$res$predict_all(n.ahead = nf, sea.on = FALSE, return.all = FALSE) # sea.on=FALSE branch too
  
  expect_equal(NROW(out0_lead$y.hat), nf)
  expect_equal(NROW(out0_targ$y.hat), nf)
  expect_equal(NROW(out0_both$y.hat), nf)
  
  # sea.period >= 2 : seasonal branches in predict_all's xpred path
  syn7_lead <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 7, xpred_logical = c(TRUE, FALSE),
                                              n_total = 80, n_est = 50)
  syn7_targ <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 7, xpred_logical = c(FALSE, TRUE),
                                              n_total = 80, n_est = 50)
  syn7_both <- make_synthetic_FilterResultsLI(resolution = "daily", sea.period = 7, xpred_logical = c(TRUE, TRUE),
                                              n_total = 80, n_est = 50)
  
  out7_lead <- syn7_lead$res$predict_all(n.ahead = nf, sea.on = TRUE, return.all = FALSE)
  out7_targ <- syn7_targ$res$predict_all(n.ahead = nf, sea.on = TRUE, return.all = FALSE)
  out7_both <- syn7_both$res$predict_all(n.ahead = nf, sea.on = TRUE, return.all = TRUE)  # return.all=TRUE branch too
  
  expect_equal(NROW(out7_lead$y.hat), nf)
  expect_equal(NROW(out7_targ$y.hat), nf)
  expect_equal(NROW(out7_both$y.hat), attr(syn7_both$res$output$model, "n") + nf)
  
  # plot_log_forecast(): xpred branch suppresses filtered level line (any(xpred_logical)==TRUE)
  expect_no_error(syn0_both$res$plot_log_forecast(Y = syn0_both$Y, n.ahead = nf))
})
