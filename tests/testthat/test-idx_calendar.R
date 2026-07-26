test_that("idx_calendar constructs with defaults and validates inputs", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days")
  expect_true(is_idx_calendar(cal))
  expect_equal(cal$pattern, 1)
  expect_equal(cal$pattern_start, 1L)
  expect_false(cal$posixct)
  expect_null(cal$anchor_name)
})

test_that("idx_calendar validates anchor_pos, amount, unit, pattern, pattern_start, posixct", {
  expect_error(idx_calendar(anchor = 0, anchor_pos = 1.5), "anchor_pos")
  expect_error(idx_calendar(anchor = 0, amount = 0), "amount")
  expect_error(idx_calendar(anchor = 0, amount = -1), "amount")
  expect_error(idx_calendar(anchor = 0, unit = 1), "unit")
  expect_error(idx_calendar(anchor = 0, pattern = c(1, -1)), "pattern")
  expect_error(idx_calendar(anchor = 0, pattern = numeric(0)), "pattern")
  expect_error(idx_calendar(anchor = 0, pattern = c(1, 1), pattern_start = 3L), "pattern_start")
  expect_error(idx_calendar(anchor = 0, posixct = NA), "posixct")
  expect_error(idx_calendar(anchor = 0, anchor_name = 5), "anchor_name")
})

test_that("is_idx_calendar distinguishes idx_calendar from other objects", {
  cal <- idx_calendar(anchor = 0)
  expect_true(is_idx_calendar(cal))
  expect_false(is_idx_calendar(list(anchor = 0)))
})

test_that("print.idx_calendar does not error", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = TRUE,
                      anchor_name = "launch")
  expect_no_error(print(cal))
})

test_that("idx_calendar_offset computes simple uniform-step offsets", {
  cal <- idx_calendar(anchor = 0, anchor_pos = 10L, amount = 1, unit = "days")
  expect_equal(idx_calendar_offset(cal, 10), 0)
  expect_equal(idx_calendar_offset(cal, 15), 5)
  expect_equal(idx_calendar_offset(cal, 5), -5)
  expect_equal(idx_calendar_offset(cal, c(5, 10, 15)), c(-5, 0, 5))
})

test_that("idx_calendar_offset respects a repeating pattern (business days)", {
  # 4 single-day steps (Mon-Thu), then a 3-day step Fri->Mon.
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days",
                      pattern = c(1, 1, 1, 1, 3), pattern_start = 1L,
                      posixct = TRUE)
  # From position 1 (Mon) to position 5 (Fri): 4 single-day steps = 4.
  expect_equal(idx_calendar_offset(cal, 5), 4)
  # From position 1 (Mon) to position 6 (next Mon): full cycle weight = 7.
  expect_equal(idx_calendar_offset(cal, 6), 7)
  # Backwards from position 1 to position -4 (previous Mon): -7.
  expect_equal(idx_calendar_offset(cal, -4), -7)
})

test_that("idx_to_date falls back to plain arithmetic for non-calendar anchors", {
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                      unit = "picoseconds")
  expect_equal(idx_to_date(cal, 1:3), c(0, 2.5, 5))
})

test_that("idx_to_date uses calendar-aware stepping for Date anchors with posixct = TRUE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = TRUE)
  out <- idx_to_date(cal, 1:5)
  expect_equal(out, seq(as.Date("2024-01-01"), by = "day", length.out = 5))
})

test_that("idx_to_date handles month/quarter/year units exactly (variable-length calendar steps)", {
  cal <- idx_calendar(anchor = as.Date("2024-01-31"), anchor_pos = 1L,
                      amount = 1, unit = "months", posixct = TRUE)
  out <- idx_to_date(cal, 1:3)
  # seq.Date with by="month" from Jan 31 correctly lands on leap-Feb 29 and Mar 31.
  expect_equal(out, seq(as.Date("2024-01-31"), by = "month", length.out = 3))
})

test_that("idx_to_date falls back to plain arithmetic when posixct is FALSE even with a Date anchor", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = FALSE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, as.Date("2024-01-01") + 0:2)
})

test_that("idx_to_date errors on non-integer number of calendar-unit steps", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 2.5, unit = "days", posixct = TRUE)
  expect_error(idx_to_date(cal, 2), "non-integer")
})

test_that("idx_to_pos is the inverse of idx_to_date for calendar-anchored daily series", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = TRUE)
  expect_equal(idx_to_pos(cal, "2024-01-10"), 10L)
  expect_equal(idx_to_pos(cal, as.Date("2024-01-01")), 1L)
})

test_that("idx_to_pos handles month/quarter/year calendar-relative units", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "quarters", posixct = TRUE)
  expect_equal(idx_to_pos(cal, as.Date("2024-10-01")), 4L)
})

test_that("idx_to_pos errors when date does not fall on a whole quarter boundary", {
  # With unit = "quarters", steps = months_diff / 3; a date exactly one
  # calendar month off the anchor's quarter cadence is not a whole number
  # of quarters away, so this trips the boundary check. (cal$amount only
  # affects the final position-spacing division, not this check, so
  # amount alone can't be used to construct a failing case here.)
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "quarters", posixct = TRUE)
  expect_error(idx_to_pos(cal, as.Date("2024-02-01")), "whole")
})

test_that("idx_to_pos treats a month/quarter/year date at calendar-month granularity (day-of-month is ignored)", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "months", posixct = TRUE)
  # Any day within the same calendar month as the anchor maps to position 1.
  expect_equal(idx_to_pos(cal, as.Date("2024-01-15")), 1L)
  expect_equal(idx_to_pos(cal, as.Date("2024-02-01")), 2L)
})

test_that("idx_to_pos errors for a non-calendar-anchored idx_calendar", {
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 1, unit = "days")
  expect_error(idx_to_pos(cal, "2024-01-10"), "not a calendar-anchored")
})

test_that("idx_offset_to_pos is the inverse of plain-arithmetic idx_to_date", {
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                      unit = "picoseconds")
  expect_equal(idx_offset_to_pos(cal, 12.5), 6L)
  expect_equal(idx_offset_to_pos(cal, 0), 1L)
})

test_that("idx_offset_to_pos errors when value does not fall on a whole step boundary", {
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                      unit = "picoseconds")
  expect_error(idx_offset_to_pos(cal, 1), "whole step boundary")
})

test_that("idx_to_date and idx_to_pos round-trip through xts_to_idx for a real dataset", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  expect_true(is_idx_series(conv$series))
  expect_true(is_idx_calendar(conv$calendar))
  
  first_date <- zoo::index(gauteng)[1]
  expect_equal(idx_to_date(conv$calendar, conv$series$start), first_date)
  expect_equal(idx_to_pos(conv$calendar, first_date), conv$series$start)
})