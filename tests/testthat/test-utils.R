test_that("df2ldl produces expected numerical outout", {
  data(england, package = 'tsgc')
  lng <- df2ldl(england$cum_cases[1:10])
  x <- as.numeric(lng[2])
  expect_equal(x,log((7915-6877)/6877))
})

test_that("df2ldl returns xts object", {
  data(england, package = 'tsgc')
  lng <- df2ldl(england$cum_cases[1:10])
  expect_s3_class(lng,"xts")
})

test_that("get_timeframe works correctly", {
  data(england, package = 'tsgc')
  x <- get_timeframe(england, '2020-03-20', '2020-03-29')
  y <- head(england,10)
  expect_equal(x,y)
})

test_that("reinitialise_dataframe works correctly", {
  data(england, package = 'tsgc')
  reinit <- reinitialise_dataframe(england,'2022-03-14')
  x <- as.numeric(reinit[1,1])
  expect_equal(x,79242)
})

test_that("Test argmax", {
  x <- xts::xts(1:10, order.by = seq(as.Date("2021-01-01"),
                                     length.out = 10, by = 1))
  x[5,] <- 20
  expect_identical(str(zoo::index(argmax(x))), str(zoo::index(x)[5]))
  expect_identical(
    str(zoo::index(argmax(x, decreasing = FALSE))),
    str(zoo::index(x)[1])
  )
})

test_that("identify error", {
  data(gauteng)
  model <- SSModelDynamicGompertz$new(Y = gauteng$cum_cases[1:100])
  res <- model$estimate()
  expect_equal(1,1)
})
