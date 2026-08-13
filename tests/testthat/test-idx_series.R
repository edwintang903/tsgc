test_that("idx_series constructs from a numeric vector", {
  x <- idx_series(c(1, 3, 6, 10), start = 5L)
  expect_true(is_idx_series(x))
  expect_equal(x$start, 5L)
  expect_equal(length(x), 4)
})

test_that("idx_series constructs from a matrix and from a data.frame", {
  m <- matrix(1:6, ncol = 2)
  xm <- idx_series(m, start = 1L)
  expect_true(is.matrix(xm$data))
  expect_equal(idx_ncol(xm), 2)
  
  df <- as.data.frame(m)
  xdf <- idx_series(df, start = 1L)
  expect_true(is.matrix(xdf$data))
  expect_equal(idx_ncol(xdf), 2)
})

test_that("idx_series rejects non-numeric data and non-integer start", {
  expect_error(idx_series(c("a", "b")), "numeric")
  expect_error(idx_series(1:5, start = 1.5), "start must be a single integer")
  expect_error(idx_series(1:5, start = c(1L, 2L)), "start must be a single integer")
})

test_that("is_idx_series distinguishes idx_series from other objects", {
  x <- idx_series(1:5)
  expect_true(is_idx_series(x))
  expect_false(is_idx_series(1:5))
  expect_false(is_idx_series(NULL))
})

test_that("as_idx_series passes through idx_series and coerces plain vectors", {
  x <- idx_series(1:5, start = 3L)
  expect_identical(as_idx_series(x), x)
  
  y <- as_idx_series(1:5, start = 2L)
  expect_true(is_idx_series(y))
  expect_equal(y$start, 2L)
  
  expect_null(as_idx_series(NULL))
})

test_that("length.idx_series returns number of observations for vector and matrix data", {
  xv <- idx_series(1:10, start = 1L)
  expect_equal(length(xv), 10)
  
  xm <- idx_series(matrix(1:20, ncol = 2), start = 1L)
  expect_equal(length(xm), 10)
})

test_that("idx_ncol and NCOL.idx_series report number of columns", {
  xv <- idx_series(1:10, start = 1L)
  expect_equal(idx_ncol(xv), 1L)
  
  xm <- idx_series(matrix(1:20, ncol = 4), start = 1L)
  expect_equal(idx_ncol(xm), 4L)
  expect_equal(NCOL.idx_series(xm), 4L)
})

test_that("idx_range and idx_positions report correct integer ranges", {
  x <- idx_series(1:10, start = 5L)
  expect_equal(idx_range(x), c(5L, 14L))
  expect_equal(idx_positions(x), 5:14)
})

test_that("[.idx_series subsets by position, not by row index", {
  x <- idx_series(101:110, start = 5L)
  sub <- x[7:9]
  expect_equal(idx_values(sub), 103:105)
  expect_equal(sub$start, 7L)
})

test_that("[.idx_series subsets matrix-valued series and supports column selection", {
  m <- matrix(1:20, ncol = 2)
  x <- idx_series(m, start = 1L)
  sub <- x[3:5]
  expect_equal(dim(idx_values(sub)), c(3, 2))
  expect_equal(sub$start, 3L)
  
  sub_col <- x[3:5, 2]
  expect_equal(as.numeric(idx_values(sub_col)), m[3:5, 2])
})

test_that("[.idx_series errors when requested positions fall outside range", {
  x <- idx_series(1:10, start = 1L)
  expect_error(x[20:22], "outside the range")
})

test_that("as.numeric/as.double/idx_values extract the raw data", {
  x <- idx_series(c(1.5, 2.5, 3.5), start = 1L)
  expect_equal(as.numeric(x), c(1.5, 2.5, 3.5))
  expect_equal(as.double(x), c(1.5, 2.5, 3.5))
  expect_equal(idx_values(x), c(1.5, 2.5, 3.5))
})

test_that("as.matrix.idx_series coerces vector and matrix data to matrix", {
  xv <- idx_series(1:5, start = 1L)
  expect_true(is.matrix(as.matrix(xv)))
  
  xm <- idx_series(matrix(1:10, ncol = 2), start = 1L)
  expect_true(is.matrix(as.matrix(xm)))
})

test_that("head.idx_series returns the first n observations by count, not by position", {
  x <- idx_series(101:110, start = 5L)
  h <- head(x, 3)
  expect_equal(idx_values(h), 101:103)
  expect_equal(h$start, 5L)
  
  expect_equal(idx_values(head(x, 100)), idx_values(x))
  expect_equal(length(head(x)), 6)
})

test_that("head.idx_series errors when n selects fewer than one observation", {
  x <- idx_series(1:10, start = 1L)
  expect_error(head(x, 0), "at least one")
})

test_that("tail.idx_series returns the last n observations by count, not by position", {
  x <- idx_series(101:110, start = 5L)
  tl <- tail(x, 3)
  expect_equal(idx_values(tl), 108:110)
  expect_equal(tl$start, 12L)
  
  expect_equal(idx_values(tail(x, 100)), idx_values(x))
  expect_equal(length(tail(x)), 6)
})

test_that("tail.idx_series errors when n selects fewer than one observation", {
  x <- idx_series(1:10, start = 1L)
  expect_error(tail(x, 0), "at least one")
})

test_that("head.idx_series and tail.idx_series work on matrix-valued series", {
  m <- matrix(1:20, ncol = 2)
  x <- idx_series(m, start = 1L)
  h <- head(x, 2)
  t <- tail(x, 2)
  expect_equal(dim(idx_values(h)), c(2, 2))
  expect_equal(dim(idx_values(t)), c(2, 2))
  expect_equal(h$start, 1L)
  expect_equal(t$start, 9L)
})

test_that("print.idx_series does not error", {
  x <- idx_series(1:5, start = 1L)
  expect_no_error(print(x))
})

test_that("idx_cbind combines series sharing the same start and length", {
  a <- idx_series(1:5, start = 1L)
  b <- idx_series(6:10, start = 1L)
  combined <- idx_cbind(a, b)
  expect_equal(idx_ncol(combined), 2)
  expect_equal(combined$start, 1L)
})

test_that("idx_cbind drops NULL arguments and returns NULL if nothing remains", {
  a <- idx_series(1:5, start = 1L)
  combined <- idx_cbind(a, NULL)
  expect_equal(idx_ncol(combined), 1)
  
  expect_null(idx_cbind(NULL, NULL))
})

test_that("idx_cbind errors when starts or lengths differ", {
  a <- idx_series(1:5, start = 1L)
  b <- idx_series(1:5, start = 2L)
  expect_error(idx_cbind(a, b), "same start and length")
  
  c <- idx_series(1:6, start = 1L)
  expect_error(idx_cbind(a, c), "same start and length")
})

test_that("idx_diff computes lagged differences and shrinks start accordingly", {
  x <- idx_series(c(1, 3, 6, 10, 15), start = 1L)
  d <- idx_diff(x, 1L)
  expect_equal(idx_values(d), c(2, 3, 4, 5))
  expect_equal(d$start, 2L)
  
  d2 <- idx_diff(x, 2L)
  expect_equal(idx_values(d2), c(5, 7, 9))
  expect_equal(d2$start, 3L)
})

test_that("idx_diff errors when lag is not smaller than series length", {
  x <- idx_series(1:5, start = 1L)
  expect_error(idx_diff(x, 5L), "lag must be smaller")
})

test_that("idx_lag shifts a series forward by k positions without altering values", {
  x <- idx_series(1:5, start = 1L)
  lagged <- idx_lag(x, 3L)
  expect_equal(idx_values(lagged), idx_values(x))
  expect_equal(lagged$start, 4L)
})

test_that("idx_rbind concatenates two adjoining series", {
  a <- idx_series(1:5, start = 1L)
  b <- idx_series(6:10, start = 6L)
  combined <- idx_rbind(a, b)
  expect_equal(idx_values(combined), 1:10)
  expect_equal(combined$start, 1L)
})

test_that("idx_rbind errors if y does not start immediately after x ends", {
  a <- idx_series(1:5, start = 1L)
  b <- idx_series(6:10, start = 7L)
  expect_error(idx_rbind(a, b), "must start immediately after")
})