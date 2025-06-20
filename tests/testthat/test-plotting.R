test_that("Test plot_new_cases() works", {
  data(gauteng, package = "tsgc")
  model <- tsgc::SSModelDynamicGompertz$new(Y = gauteng[1:50], q = 0.005)
  res <- model$estimate()
  tsgc::plot_new_cases(res)
  expect_equal(1,1)
})

test_that("Test plot_log_forecast() works", {
  data(gauteng, package = "tsgc")
  model <- tsgc::SSModelDynamicGompertz$new(Y = gauteng[1:50], q = 0.005)
  res <- model$estimate()
  tsgc::plot_log_forecast(res, gauteng)
  expect_equal(1,1)
})

test_that("Test plot_gy_components() works", {
  data(gauteng, package = "tsgc")
  model <- tsgc::SSModelDynamicGompertz$new(Y = gauteng[1:50], q = 0.005)
  res <- model$estimate()
  tsgc::plot_gy_components(res)
  expect_equal(1,1)
})


test_that("Test plot_gy_ci() works", {
  data(gauteng, package = "tsgc")
  model <- tsgc::SSModelDynamicGompertz$new(Y = gauteng[1:50], q = 0.005)
  res <- model$estimate()
  tsgc::plot_gy_ci(res)
  expect_equal(1,1)
})


test_that("Test plot_holdout() works", {
  data(gauteng, package = "tsgc")
  model <- tsgc::SSModelDynamicGompertz$new(Y = gauteng[1:50], q = 0.005)
  idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
  model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
  res <- model$estimate()
  # Plot forecasts and outcomes over evaluation period
  plot_holdout(res = res, Y = gauteng)
  expect_equal(1,1)
})
