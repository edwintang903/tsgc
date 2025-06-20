test_that("Test Gauteng data is correct class and type", {
  data("gauteng")
  expect_type(gauteng, 'integer')
  expect_s3_class(gauteng, 'xts')
})

test_that("Test England data is correct class", {
  data("england")
  expect_s3_class(england, 'xts')
})

test_that("Test UK and Italy data is correct class and type", {
  data("ukitaly")
  expect_type(ukitaly, 'integer')
  expect_s3_class(ukitaly, 'xts')
})

test_that("Test Gauteng weather data is correct class", {
  data("gauteng_weather_2021")
  expect_s3_class(gauteng_weather_2021, 'xts')
})

test_that("Test England weather data is correct class", {
  data("england_weather_2021")
  expect_s3_class(england_weather_2021, 'xts')
})
