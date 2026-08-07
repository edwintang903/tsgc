test_that("xts_to_idx converts a daily xts object to an idx_series + idx_calendar", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  
  expect_true(is_idx_series(conv$series))
  expect_true(is_idx_calendar(conv$calendar))
  expect_equal(as.numeric(idx_values(conv$series)), as.numeric(zoo::coredata(england$cum_cases)))
  expect_equal(conv$calendar$anchor, zoo::index(england$cum_cases)[1])
})

test_that("xts_to_idx honours a custom start.pos for alignment", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases, start.pos = 100L)
  expect_equal(conv$series$start, 100L)
  expect_equal(conv$calendar$anchor_pos, 100L)
})

test_that("xts_to_idx detects a plain daily pattern", {
  x <- xts::xts(cumsum(rpois(30, 5)) + 1, order.by = as.Date("2024-01-01") + 0:29)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$amount, 1)
  expect_equal(conv$calendar$unit, "days")
  expect_equal(conv$calendar$pattern, 1)
  expect_true(conv$calendar$posixct)
})

test_that("xts_to_idx detects a business-day pattern (weekends skipped)", {
  # 41 business days from Mon 1 Jan gives 40 gaps, a whole multiple of
  # the 5-day Mon-Fri cycle length that find_pattern() requires.
  bdays <- seq(as.Date("2024-01-01"), by = "day", length.out = 60)
  bdays <- bdays[!weekdays(bdays) %in% c("Saturday", "Sunday")][1:41]
  x <- xts::xts(cumsum(rpois(length(bdays), 5)) + 1, order.by = bdays)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$unit, "days")
  expect_equal(conv$calendar$pattern, c(1, 1, 1, 1, 3))
})

test_that("xts_to_idx detects a monthly Date index as amount=1, unit='months'", {
  months <- seq(as.Date("2024-01-01"), by = "month", length.out = 24)
  x <- xts::xts(cumsum(rpois(24, 5)) + 1, order.by = months)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$amount, 1)
  expect_equal(conv$calendar$unit, "months")
})

test_that("xts_to_idx detects a quarterly Date index as amount=1, unit='quarters'", {
  qtrs <- seq(as.Date("2024-01-01"), by = "quarter", length.out = 12)
  x <- xts::xts(cumsum(rpois(12, 5)) + 1, order.by = qtrs)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$unit, "quarters")
})

test_that("xts_to_idx detects a yearly Date index as amount=1, unit='years'", {
  yrs <- seq(as.Date("2024-01-01"), by = "year", length.out = 8)
  x <- xts::xts(cumsum(rpois(8, 5)) + 1, order.by = yrs)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$unit, "years")
})

test_that("xts_to_idx detects a yearqtr index", {
  qtrs <- zoo::as.yearqtr("2020 Q1") + 0:9 * 0.25
  x <- xts::xts(cumsum(rpois(10, 5)) + 1, order.by = qtrs)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$unit, "quarters")
  expect_true(inherits(conv$calendar$anchor, "yearqtr"))
})

test_that("xts_to_idx detects a yearmon index", {
  mons <- zoo::as.yearmon("2020-01") + 0:11 * (1/12)
  x <- xts::xts(cumsum(rpois(12, 5)) + 1, order.by = mons)
  conv <- xts_to_idx(x)
  expect_equal(conv$calendar$unit, "months")
  expect_true(inherits(conv$calendar$anchor, "yearmon"))
})

test_that("xts_to_idx errors on an irregular index when detect = TRUE", {
  gaps <- c(1, 2, 1, 6, 1, 3, 2, 5, 1, 1, 4, 2)
  irregular_dates <- as.Date("2024-01-01") + c(0, cumsum(gaps))
  x <- xts::xts(seq_along(irregular_dates), order.by = irregular_dates)
  expect_error(xts_to_idx(x), "could not detect")
})

test_that("xts_to_idx errors on duplicate index values", {
  d <- as.Date("2024-01-01") + c(0, 1, 1, 2)
  x <- xts::xts(1:4, order.by = d)
  expect_error(xts_to_idx(x), "duplicate")
})

test_that("xts_to_idx sorts an out-of-order index before conversion", {
  d <- as.Date("2024-01-01") + c(2, 0, 1)
  x <- xts::xts(c(30, 10, 20), order.by = d)
  conv <- xts_to_idx(x)
  expect_equal(as.numeric(idx_values(conv$series)), c(10, 20, 30))
  expect_equal(conv$calendar$anchor, as.Date("2024-01-01"))
})

test_that("xts_to_idx with detect = FALSE uses the supplied amount/unit rather than detecting", {
  x <- xts::xts(1:5, order.by = as.Date("2024-01-01") + 0:4)
  conv <- xts_to_idx(x, detect = FALSE, amount = 3, unit = "hours")
  expect_equal(conv$calendar$amount, 3)
  expect_equal(conv$calendar$unit, "hours")
  expect_equal(conv$calendar$pattern, 1)
})

test_that("xts_to_idx detect = TRUE errors with fewer than 2 rows", {
  x <- xts::xts(1, order.by = as.Date("2024-01-01"))
  expect_error(xts_to_idx(x), "at least 2 rows")
})

test_that("df2ldl produces expected numerical output", {
  data(england, package = "tsgc")
  conv <- xts_to_idx(england$cum_cases)
  x10 <- get_timeframe(conv$series, conv$series$start, conv$series$start + 9)
  lng <- df2ldl(x10)
  expect_equal(idx_values(lng)[2], log((7915 - 6877) / 6877))
})

test_that("df2ldl returns an idx_series with a leading NA and errors on bad input", {
  x <- idx_series(cumsum(1:10) + 1, start = 1L)
  lng <- df2ldl(x)
  expect_true(is_idx_series(lng))
  expect_true(is.na(idx_values(lng)[1]))
  expect_equal(length(lng), length(x))
  
  expect_error(df2ldl(1:10), "idx_series")
  expect_error(df2ldl(idx_series(matrix(1:10, ncol = 2))), "1 data column")
})

test_that("df2ldl errors on negative levels or genuinely negative increments", {
  x_neg <- idx_series(c(10, -5, 20, 30), start = 1L)
  expect_error(df2ldl(x_neg), "negative")
  
  x_dec <- idx_series(c(10, 20, 15, 30), start = 1L)
  expect_error(df2ldl(x_dec), "nonpositive increments")
})

test_that("df2ldl does not error on an exact plateau (zero increment), but produces -Inf at that position", {
  x_flat <- idx_series(c(10, 20, 20, 30), start = 1L)
  ldl <- df2ldl(x_flat)
  expect_true(is.infinite(idx_values(ldl)[3]) && idx_values(ldl)[3] < 0)
})

test_that("get_timeframe subsets an idx_series by integer position range", {
  x <- idx_series(1:30, start = 1L)
  sub <- get_timeframe(x, 5, 10)
  expect_equal(idx_values(sub), 5:10)
  expect_equal(sub$start, 5L)
})

test_that("get_timeframe defaults end to the last position in df", {
  x <- idx_series(1:30, start = 1L)
  sub <- get_timeframe(x, 25)
  expect_equal(idx_values(sub), 25:30)
})

test_that("get_timeframe clamps start/end to the available range of df", {
  x <- idx_series(1:10, start = 1L)
  sub <- get_timeframe(x, -5, 20)
  expect_equal(idx_values(sub), 1:10)
})

test_that("get_timeframe errors when start is after end within the available range", {
  x <- idx_series(1:10, start = 1L)
  expect_error(get_timeframe(x, 9, 3), "start is after end")
})

test_that("get_timeframe returns NULL when df is NULL, and errors on non-idx_series input", {
  expect_null(get_timeframe(NULL, 1, 10))
  expect_error(get_timeframe(1:10, 1, 10), "idx_series")
})

test_that("reinitialise_dataframe rebases a series relative to the value just before reinit.idx", {
  x <- idx_series(c(100, 110, 125, 140, 160), start = 1L)
  reinit <- reinitialise_dataframe(x, 3L)
  expect_equal(reinit$start, 3L)
  expect_equal(idx_values(reinit), idx_values(x)[3:5] - idx_values(x)[2])
  expect_equal(idx_values(reinit)[1], 15)
})

test_that("reinitialise_dataframe errors when reinit.idx is out of range or has no preceding value", {
  x <- idx_series(1:10, start = 1L)
  expect_error(reinitialise_dataframe(x, 1L), "not present")
  expect_error(reinitialise_dataframe(x, 15L), "not present")
})

test_that("argmax returns expected output for an idx_series", {
  x <- idx_series(c(1, 2, 20, 3, 4), start = 1L)
  mx <- argmax(x)
  expect_true(is_idx_series(mx))
  expect_equal(idx_values(mx), 20)
  expect_equal(mx$start, 3L)
  
  mn <- argmax(x, decreasing = FALSE)
  expect_equal(idx_values(mn), 1)
  expect_equal(mn$start, 1L)
})

test_that("argmax returns expected output for a plain numeric vector", {
  x <- c(5, 1, 9, 3)
  expect_equal(argmax(x), 9)
  expect_equal(argmax(x, decreasing = FALSE), 1)
})

test_that("mapes works and returns five named error metrics", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(
    Y = conv$series, q = 0.005,
    start = est.start, end = est.end
  )
  res <- estimate(model)
  
  error_metrics <- mapes(res, n.ahead = 7, Y = conv$series)
  
  expect_type(error_metrics, "list")
  expected_names <- c("mape", "smape", "mae", "rmse", "coverage")
  expect_equal(sort(names(error_metrics)), sort(expected_names))
})

test_that("estimate_r0 returns positions and calendar dates in a data frame", {
  set.seed(1)
  Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1L)
  cal <- idx_calendar(as.Date("2021-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days")
  model <- SSModelDynamicGompertz$new(
    Y = Y, q = NULL, end = 100, calendar = cal
  )
  res <- estimate(model)
  
  r_t <- estimate_r0(res, gen_int = 5, n.ahead = 7)
  
  expect_s3_class(r_t, "data.frame")
  expect_equal(names(r_t), c("Date", "Position", "fit", "lower", "upper"))
  expect_equal(nrow(r_t), 7)
  expect_equal(r_t$Date, idx_to_date(cal, r_t$Position))
  expect_true(all(r_t$lower <= r_t$fit))
  expect_true(all(r_t$fit <= r_t$upper))
})

test_that("estimate_r0 retains positions when no calendar is available", {
  set.seed(2)
  Y <- idx_series(cumsum(rpois(80, 8)) + 1, start = 1L)
  res <- estimate(SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 60))
  
  r_t <- estimate_r0(res, gen_int = 5, n.ahead = 5)
  
  expect_equal(names(r_t), c("Position", "fit", "lower", "upper"))
  expect_equal(nrow(r_t), 5)
})

test_that("estimate_r0 errors when res is not a FilterResults or FilterResultsLI object", {
  expect_error(estimate_r0(list(), gen_int = 5), "FilterResults")
})

test_that("add_daily_ldl handles LeadIndCol logic correctly", {
  # Column 1: 10, 20, 40, 70 (increments: NA, 10, 20, 30)
  # Column 2: 5, 15, 25, 35  (increments: NA, 10, 10, 10)
  data_mat <- cbind(col_a = c(10, 20, 40, 70), col_b = c(5, 15, 25, 35))
  x <- idx_series(data_mat, start = 1L)
  
  res1 <- add_daily_ldl(x, LeadIndCol = 1)
  expect_equal(idx_values(res1$cLead), c(10, 20, 40, 70))
  expect_equal(idx_values(res1$cTarg), c(5, 15, 25, 35))
  
  res2 <- add_daily_ldl(x, LeadIndCol = 2)
  expect_equal(idx_values(res2$cLead), c(5, 15, 25, 35))
  expect_equal(idx_values(res2$cTarg), c(10, 20, 40, 70))
})

test_that("add_daily_ldl computes increments correctly", {
  data_mat <- cbind(col_a = c(10, 20, 40, 70), col_b = c(5, 15, 25, 35))
  x <- idx_series(data_mat, start = 1L)
  
  res <- add_daily_ldl(x, LeadIndCol = 1)
  
  expect_equal(idx_values(res$newLead), c(10, 20, 30))
  expect_equal(idx_values(res$newTarg), c(10, 10, 10))
  expect_equal(res$newLead$start, 2L)
})

test_that("add_daily_ldl output structure is correct", {
  data_mat <- cbind(col_a = c(10, 20, 40, 70), col_b = c(5, 15, 25, 35))
  x <- idx_series(data_mat, start = 1L)
  
  res <- add_daily_ldl(x)
  
  expect_equal(names(res), c("cLead", "cTarg", "newLead", "newTarg", "LDLlead", "LDLtarg"))
  expect_true(all(vapply(res, is_idx_series, logical(1))))
})

test_that("add_daily_ldl errors on non-idx_series or wrong number of columns", {
  expect_error(add_daily_ldl(1:10), "idx_series")
  expect_error(add_daily_ldl(idx_series(1:10)), "exactly two series")
  expect_error(add_daily_ldl(idx_series(matrix(1:10, ncol = 2)), LeadIndCol = 3), "LeadIndCol")
})

test_that("cross_val reports the correct criterion (e.g., 'rmse')", {
  data(ukitaly, package = "tsgc")
  conv <- xts_to_idx(ukitaly)
  Yuk <- idx_series(idx_values(conv$series)[, "UK"], start = conv$series$start)
  
  est.start <- idx_to_pos(conv$calendar, as.Date("2020-02-25"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2020-04-01"))
  
  cv_models <- list()
  cv_models[["Vanilla_q"]] <- SSModelDynamicGompertz$new(Y = Yuk, q = 0.005, start = est.start, end = est.end)
  cv_models[["Vanilla_ar1"]] <- SSModelDynamicGompertz$new(Y = Yuk, start = est.start, end = est.end, ar1 = TRUE)
  cv_models[["Lag7"]] <- SSModelLeadingIndicator$new(Y = conv$series, start = est.start, end = est.end, n.lag = 7)
  
  n.ahead <- 7
  
  cv_result_rmse <- cross_val(
    Y = conv$series,
    model_list = cv_models,
    est.end = est.end,
    n.ahead = n.ahead,
    n.estimate = 1,
    gap = 2,
    criterion = "rmse"
  )
  
  col_name <- as.character(est.end)
  expect_type(cv_result_rmse[[col_name]], "double")
  expect_true(all(cv_result_rmse[[col_name]] > 0))
  
  model_q_test <- SSModelDynamicGompertz$new(Y = Yuk, q = 0.005, start = est.start, end = est.end)
  res_q_test <- estimate(model_q_test)
  expected_rmse <- round(mapes(res_q_test, n.ahead, Yuk)[["rmse"]], 2)
  
  expect_equal(cv_result_rmse[cv_result_rmse$Model == "Vanilla_q", col_name][[1]], expected_rmse)
})

test_that("cross_val errors on non-idx_series Y, bad est.end, or non-positive n.ahead", {
  x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1L)
  model_list <- list(m1 = SSModelDynamicGompertz$new(Y = x, end = 20))
  
  expect_error(cross_val(Y = 1:10, model_list = model_list, est.end = 20), "idx_series")
  expect_error(cross_val(Y = x, model_list = model_list, est.end = c(1, 2)), "single integer")
  expect_error(cross_val(Y = x, model_list = model_list, est.end = 20, n.ahead = 0), "positive")
})

test_that("write_results writes forecast, filtered level/slope and growth rate CSVs to disk", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, start = est.start, end = est.end)
  res <- estimate(model)
  
  res.dir <- tempdir()
  prefix <- paste0("test_", as.integer(Sys.time()), "_")
  
  expect_no_error(
    write_results(res = res, res.dir = res.dir, n.ahead = 5, prefix = prefix, confidence.level = 0.68)
  )
  
  expect_true(file.exists(file.path(res.dir, paste0(prefix, "cases_fcst.csv"))))
  expect_true(file.exists(file.path(res.dir, paste0(prefix, "trend_slope_filt.csv"))))
  expect_true(file.exists(file.path(res.dir, paste0(prefix, "log_gr_level_filt.csv"))))
  expect_true(file.exists(file.path(res.dir, paste0(prefix, "cases_gr.csv"))))
})

test_that("write_results errors when res is not a FilterResults or FilterResultsLI object", {
  expect_error(
    write_results(res = list(), res.dir = tempdir(), n.ahead = 5),
    "FilterResults"
  )
})

test_that("check_variance_boundary flags a scalar variance at or below the tolerance", {
  expect_true(check_variance_boundary(0, boundary_tol = 1e-6))
  expect_true(check_variance_boundary(1e-7, boundary_tol = 1e-6))
  expect_false(check_variance_boundary(0.01, boundary_tol = 1e-6))
})

test_that("check_variance_boundary works on a matrix-valued variance parameter (e.g. leading-indicator H) without erroring, unlike a direct && on a multi-element value", {
  H_scalar <- matrix(0.05, 1, 1)
  H_matrix_ok <- matrix(c(0.05, 0.001, 0.001, 0.03), nrow = 2)
  H_matrix_boundary <- matrix(c(0.05, 0.001, 0.001, 1e-8), nrow = 2)
  
  expect_false(check_variance_boundary(H_scalar, boundary_tol = 1e-6))
  expect_false(check_variance_boundary(H_matrix_ok, boundary_tol = 1e-6))
  expect_true(check_variance_boundary(H_matrix_boundary, boundary_tol = 1e-6))
})

test_that("check_variance_boundary treats NULL/NA arguments as absent rather than erroring or matching", {
  expect_false(check_variance_boundary(NULL, boundary_tol = 1e-6))
  expect_false(check_variance_boundary(NA, boundary_tol = 1e-6))
  expect_true(check_variance_boundary(NULL, 1e-8, boundary_tol = 1e-6))
})

test_that("print_model_diagnostics runs without error on a FilterResults object and reports log-likelihood and a boundary check", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng$cum_cases)
  est.start <- idx_to_pos(conv$calendar, as.Date("2021-02-01"))
  est.end   <- idx_to_pos(conv$calendar, as.Date("2021-04-19"))
  
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, start = est.start, end = est.end)
  res <- estimate(model)
  
  out <- capture.output(print_model_diagnostics(res))
  expect_true(any(grepl("Log-likelihood", out)))
  expect_true(any(grepl("Boundary check", out)))
  expect_true(any(grepl("Recursive residuals", out)))
})

test_that("print_model_diagnostics runs without error on a FilterResultsLI object, where H is a 2x2 matrix rather than a scalar", {
  lead <- seq(100, 400, length.out = 40)
  targ <- seq(50, 200, length.out = 40)
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 4, sea.period = 0, LeadIndCol = 1)
  res <- mod$estimate()
  
  expect_no_error(out <- capture.output(print_model_diagnostics(res)))
  expect_true(any(grepl("Log-likelihood", out)))
  expect_true(any(grepl("Boundary check", out)))
})

test_that("print_model_diagnostics errors clearly when res is not a FilterResults or FilterResultsLI object", {
  expect_error(print_model_diagnostics(list()), "FilterResults")
})

test_that("FilterResultsLI's summary() prints the 2x2 observation-noise matrix as a labelled matrix, not flattened into a single cat() line", {
  lead <- seq(100, 400, length.out = 40)
  targ <- seq(50, 200, length.out = 40)
  Y_li <- idx_series(cbind(lead_col = lead, targ_col = targ), start = 1L)
  
  mod <- SSModelLeadingIndicator$new(Y = Y_li, n.lag = 4, sea.period = 0, LeadIndCol = 1)
  res <- mod$estimate()
  
  out <- capture.output(res$summary())
  expect_true(any(grepl("Observation equation noise", out)))
  noise_line <- grep("Observation equation noise", out)
  expect_true(any(grepl("\\[1,\\]", out[(noise_line + 1):length(out)][1:3])))
})