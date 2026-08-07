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
  expect_equal(out, seq(as.Date("2024-01-31"), by = "month", length.out = 3))
})

test_that("idx_to_date supports a yearqtr anchor with unit = 'quarters'", {
  cal <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                      amount = 1, unit = "quarters", posixct = TRUE)
  out <- idx_to_date(cal, 1:4)
  expect_equal(out, zoo::as.yearqtr("2024 Q1") + 0:3 * 0.25)
  expect_true(inherits(out, "yearqtr"))
})

test_that("idx_to_date supports a yearqtr anchor with unit = 'years'", {
  cal <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                      amount = 1, unit = "years", posixct = TRUE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, zoo::as.yearqtr("2024 Q1") + 0:2)
})

test_that("idx_to_date supports a yearmon anchor with unit = 'months'", {
  cal <- idx_calendar(anchor = zoo::as.yearmon("2024-01"), anchor_pos = 1L,
                      amount = 1, unit = "months", posixct = TRUE)
  out <- idx_to_date(cal, 1:5)
  expect_equal(out, zoo::as.yearmon("2024-01") + 0:4 * (1/12))
  expect_true(inherits(out, "yearmon"))
})

test_that("idx_to_date supports a yearmon anchor with unit = 'years'", {
  cal <- idx_calendar(anchor = zoo::as.yearmon("2024-01"), anchor_pos = 1L,
                      amount = 1, unit = "years", posixct = TRUE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, zoo::as.yearmon("2024-01") + 0:2)
})

test_that("idx_to_date respects amount for yearqtr/yearmon anchors (e.g. every 2 quarters)", {
  cal <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                      amount = 2, unit = "quarters", posixct = TRUE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, zoo::as.yearqtr("2024 Q1") + c(0, 2, 4) * 0.25)
})

test_that("idx_to_date rejects a 'months' unit for a yearqtr anchor", {
  cal <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                      amount = 1, unit = "months", posixct = TRUE)
  expect_error(idx_to_date(cal, 1:3), "not supported")
})

test_that("idx_to_date rejects a 'quarters' unit for a yearmon anchor", {
  cal <- idx_calendar(anchor = zoo::as.yearmon("2024-01"), anchor_pos = 1L,
                      amount = 1, unit = "quarters", posixct = TRUE)
  expect_error(idx_to_date(cal, 1:3), "not supported")
})

test_that("idx_to_date rejects a sub-resolution unit (days) for yearqtr/yearmon anchors", {
  cal_qtr <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                          amount = 1, unit = "days", posixct = TRUE)
  cal_mon <- idx_calendar(anchor = zoo::as.yearmon("2024-01"), anchor_pos = 1L,
                          amount = 1, unit = "days", posixct = TRUE)
  expect_error(idx_to_date(cal_qtr, 1:3), "not supported")
  expect_error(idx_to_date(cal_mon, 1:3), "not supported")
})

test_that("idx_to_date falls back to plain arithmetic for yearqtr/yearmon anchors when posixct = FALSE", {
  cal_qtr <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                          amount = 1, unit = "quarters", posixct = FALSE)
  cal_mon <- idx_calendar(anchor = zoo::as.yearmon("2024-01"), anchor_pos = 1L,
                          amount = 1, unit = "months", posixct = FALSE)
  # yearqtr/yearmon are numeric under the hood, so plain arithmetic still
  # lands close to the right value here, via the generic branch.
  expect_equal(idx_to_date(cal_qtr, 1:3), zoo::as.yearqtr("2024 Q1") + 0:2)
  expect_equal(idx_to_date(cal_mon, 1:3), zoo::as.yearmon("2024-01") + 0:2)
})

test_that("idx_to_date respects a repeating pattern for a yearqtr anchor", {
  # 3rd slot (Q3) has weight 2: stepping onto/past it advances 2 quarters.
  cal <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                      amount = 1, unit = "quarters",
                      pattern = c(1, 1, 2, 1), pattern_start = 1L,
                      posixct = TRUE)
  out <- idx_to_date(cal, 1:4)
  # Pattern-weighted offsets at positions 1:4 are 0, 1, 2, 4.
  expect_equal(out, zoo::as.yearqtr("2024 Q1") + c(0, 1, 2, 4) * 0.25)
})

test_that("idx_to_date falls back to plain arithmetic when posixct is FALSE even with a Date anchor", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = FALSE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, as.Date("2024-01-01") + 0:2)
})

test_that("idx_to_date errors on non-integer number of calendar-unit steps (fractional pattern weights)", {
  # A fractional pattern can produce a non-integer offset from an
  # otherwise whole integer position.
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days",
                      pattern = c(0.5, 1.5), pattern_start = 1L,
                      posixct = TRUE)
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
  # A date one calendar month off the anchor's quarter cadence is not a
  # whole number of quarters away.
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

test_that("idx_to_pos is the inverse of idx_to_date for a yearqtr anchor", {
  cal <- idx_calendar(anchor = zoo::as.yearqtr("2024 Q1"), anchor_pos = 1L,
                      amount = 1, unit = "quarters", posixct = TRUE)
  expect_equal(idx_to_pos(cal, zoo::as.yearqtr("2024 Q4")), 4L)
  expect_equal(idx_to_pos(cal, "2024 Q1"), 1L)
})

test_that("idx_to_pos is the inverse of idx_to_date for a yearmon anchor", {
  cal <- idx_calendar(anchor = zoo::as.yearmon("2024-01"), anchor_pos = 1L,
                      amount = 1, unit = "months", posixct = TRUE)
  expect_equal(idx_to_pos(cal, zoo::as.yearmon("2024-05")), 5L)
})

test_that("idx_to_pos: (via idx_step_date_to_offset) binary-searches the inverse of an idx_calendar_step calendar", {
  cal <- idx_calendar_step(
    anchor = as.POSIXct("2024-01-01 00:00:03", tz = "UTC"),
    step = idx_step(quarters = 1, seconds = 3), anchor_pos = 1L, posixct = TRUE
  )
  expect_equal(idx_to_pos(cal, as.POSIXct("2024-01-01 00:00:03", tz = "UTC")), 1L)
  expect_equal(idx_to_pos(cal, as.POSIXct("2024-10-01 00:00:12", tz = "UTC")), 4L)
})

test_that("idx_to_pos: idx_step_date_to_offset errors when date is not on a step boundary", {
  cal <- idx_calendar_step(anchor = as.Date("2024-01-01"),
                           step = idx_step(months = 1), anchor_pos = 1L, posixct = TRUE)
  expect_error(idx_to_pos(cal, as.Date("2024-01-15")), "whole step boundary")
})

test_that("idx_to_pos: idx_step_date_to_offset handles negative offsets (date before anchor)", {
  cal <- idx_calendar_step(anchor = as.Date("2024-06-01"),
                           step = idx_step(months = 1), anchor_pos = 10L, posixct = TRUE)
  expect_equal(idx_to_pos(cal, as.Date("2024-01-01")), 5L)
})

test_that("idx_to_pos: (via idx_multi_step_date_to_offset) walks the inverse of a multi_step_pattern calendar forward and backward", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                                 multi_step = msp, posixct = TRUE)
  expect_equal(idx_to_pos(cal, as.Date("2024-01-07")), 3L)
  expect_equal(idx_to_pos(cal, as.Date("2024-02-10")), 5L)
  expect_equal(idx_to_pos(cal, as.Date("2024-01-01")), 1L)
})

test_that("idx_to_pos: idx_multi_step_date_to_offset errors when date does not fall on a step boundary", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                                 multi_step = msp, posixct = TRUE)
  expect_error(idx_to_pos(cal, as.Date("2024-01-02")), "step boundary")
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

test_that("idx_step constructs with defaults and validates fields", {
  s <- idx_step(days = 1)
  expect_true(is_idx_step(s))
  expect_equal(s$days, 1)
  expect_equal(s$years, 0)
  expect_error(idx_step(), "non-zero")
  expect_error(idx_step(days = -1), "days")
  expect_error(idx_step(seconds = c(1, 2)), "seconds")
})

test_that("idx_step_needs_posixct is TRUE only when hours/minutes/seconds are set", {
  expect_false(idx_step_needs_posixct(idx_step(months = 1)))
  expect_true(idx_step_needs_posixct(idx_step(seconds = 3)))
  expect_true(idx_step_needs_posixct(idx_step(quarters = 1, seconds = 3)))
})

test_that("print.idx_step does not error and omits zero fields", {
  s <- idx_step(months = 1, seconds = 13)
  expect_no_error(print(s))
})

test_that("idx_step_add applies a compound step (quarters + seconds) to a POSIXct anchor", {
  anchor <- as.POSIXct("2024-01-01 00:00:03", tz = "UTC")
  s <- idx_step(quarters = 1, seconds = 3)
  out <- idx_step_add(anchor, s, 0:3)
  # Each field is scaled by n independently: n quarters (calendar-aware)
  # plus n * 3 seconds (cumulative, not a one-time nudge).
  expect_equal(out[1], anchor)
  expect_equal(out[2], as.POSIXct("2024-04-01 00:00:06", tz = "UTC"))
  expect_equal(out[3], as.POSIXct("2024-07-01 00:00:09", tz = "UTC"))
  expect_equal(out[4], as.POSIXct("2024-10-01 00:00:12", tz = "UTC"))
})

test_that("idx_step_add applies a compound step (months + seconds) to a POSIXct anchor", {
  anchor <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  s <- idx_step(months = 1, seconds = 13)
  out <- idx_step_add(anchor, s, 0:2)
  expect_equal(out[1], anchor)
  expect_equal(out[2], as.POSIXct("2024-02-01 00:00:13", tz = "UTC"))
  expect_equal(out[3], as.POSIXct("2024-03-01 00:00:26", tz = "UTC"))
})

test_that("idx_step_add keeps components independent rather than normalizing (months + days)", {
  anchor <- as.Date("2024-01-01")
  s <- idx_step(months = 1, days = 31)
  out <- idx_step_add(anchor, s, 1)
  # "One month, then 31 days" - not simplified/normalized to two months.
  expect_equal(out, as.Date("2024-02-01") + 31)
})

test_that("idx_step_add errors when a sub-day component is used with a yearqtr/yearmon anchor", {
  s <- idx_step(quarters = 1, seconds = 3)
  expect_error(idx_step_add(zoo::as.yearqtr("2024 Q1"), s, 1), "yearqtr/yearmon")
})

test_that("idx_step_add errors for a 'months' component on a yearqtr anchor, and vice versa", {
  expect_error(idx_step_add(zoo::as.yearqtr("2024 Q1"), idx_step(months = 1), 1), "yearqtr")
  expect_error(idx_step_add(zoo::as.yearmon("2024-01"), idx_step(quarters = 1), 1), "yearmon")
})

test_that("idx_step_add handles negative n (stepping backward)", {
  anchor <- as.Date("2024-03-01")
  s <- idx_step(months = 1)
  out <- idx_step_add(anchor, s, -2)
  expect_equal(out, as.Date("2024-01-01"))
})

test_that("idx_calendar_step constructs a calendar with a compound step and NA amount/unit", {
  cal <- idx_calendar_step(
    anchor = as.POSIXct("2024-01-01 00:00:03", tz = "UTC"),
    step = idx_step(quarters = 1, seconds = 3), posixct = TRUE
  )
  expect_true(is_idx_calendar(cal))
  expect_true(is.na(cal$amount))
  expect_true(is_idx_step(cal$step))
})

test_that("idx_calendar_step errors when step is not an idx_step", {
  expect_error(idx_calendar_step(anchor = as.Date("2024-01-01"), step = "not a step"), "idx_step")
})

test_that("idx_to_date delegates to idx_step_add for an idx_calendar_step calendar", {
  cal <- idx_calendar_step(
    anchor = as.POSIXct("2024-01-01 00:00:03", tz = "UTC"),
    step = idx_step(quarters = 1, seconds = 3), posixct = TRUE
  )
  out <- idx_to_date(cal, 1:4)
  expect_equal(out[1], as.POSIXct("2024-01-01 00:00:03", tz = "UTC"))
  expect_equal(out[4], as.POSIXct("2024-10-01 00:00:12", tz = "UTC"))
})

test_that("idx_to_date errors for an idx_calendar_step calendar with posixct = FALSE (no single amount to fall back to)", {
  cal <- idx_calendar_step(anchor = as.Date("2024-01-01"),
                           step = idx_step(months = 1, days = 1), posixct = FALSE)
  expect_error(idx_to_date(cal, 1:3), "no single 'amount'")
})

test_that("idx_calendar's primary constructor still builds an equivalent internal idx_step", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "months", posixct = TRUE)
  expect_true(is_idx_step(cal$step))
  expect_equal(cal$step$months, 1)
  expect_null(cal$multi_step)
})

test_that("multi_step_pattern constructs from idx_steps and validates its arguments", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  expect_true(is_multi_step_pattern(msp))
  expect_equal(length(msp), 3)
  expect_error(multi_step_pattern(), "at least one")
  expect_error(multi_step_pattern(idx_step(days = 1), "not a step"), "idx_step objects")
})

test_that("print.multi_step_pattern does not error", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(months = 1))
  expect_no_error(print(msp))
})

test_that("idx_calendar_multi_step constructs and validates pattern_start", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), multi_step = msp, posixct = TRUE)
  expect_true(is_idx_calendar(cal))
  expect_true(is_multi_step_pattern(cal$multi_step))
  expect_error(
    idx_calendar_multi_step(anchor = as.Date("2024-01-01"), multi_step = msp,
                            pattern_start = 4L, posixct = TRUE),
    "pattern_start"
  )
  expect_error(
    idx_calendar_multi_step(anchor = as.Date("2024-01-01"), multi_step = "not a pattern"),
    "multi_step_pattern"
  )
})

test_that("idx_to_date walks a multi_step_pattern (3 days, 3 days, 1 month) forward correctly", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                                 multi_step = msp, posixct = TRUE)
  out <- idx_to_date(cal, 1:5)
  expect_equal(out, as.POSIXct(c("2024-01-01", "2024-01-04", "2024-01-07",
                                 "2024-02-07", "2024-02-10"), tz = attr(out, "tzone")))
})

test_that("idx_to_date walks a multi_step_pattern backward from a later anchor", {
  # Walking backward from position 5 undoes the slot before the current
  # pointer, so it doesn't retrace the forward walk's slot sequence;
  # here it walks slots 3, 2, 1, 3 undone in turn.
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-02-10"), anchor_pos = 5L,
                                 multi_step = msp, posixct = TRUE)
  out <- idx_to_date(cal, 1:5)
  expect_equal(out, as.POSIXct(c("2023-12-04", "2024-01-04", "2024-01-07",
                                 "2024-01-10", "2024-02-10"), tz = attr(out, "tzone")))
})

test_that("idx_to_date walks a multi_step_pattern backward as the exact inverse of the forward walk, when pattern_start is set to match", {
  # pattern_start must identify the slot whose step leads OUT of
  # anchor_pos in the forward direction for the walks to match exactly.
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  cal_fwd <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                                     multi_step = msp, posixct = TRUE)
  fwd <- idx_to_date(cal_fwd, 1:5)
  
  cal_bwd <- idx_calendar_multi_step(anchor = as.Date("2024-02-10"), anchor_pos = 5L,
                                     multi_step = msp, pattern_start = 2L, posixct = TRUE)
  out <- idx_to_date(cal_bwd, 1:5)
  expect_equal(out, fwd)
})

test_that("idx_to_date respects pattern_start for a multi_step_pattern calendar", {
  # Same 3-slot cycle, but anchor_pos sits at the 2nd slot (i.e. the step
  # out of anchor is the second 'days=3' step, not the first).
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                                 multi_step = msp, pattern_start = 2L, posixct = TRUE)
  out <- idx_to_date(cal, 1:4)
  # pos1 = anchor; pos2 = +3 days (slot 2); pos3 = +1 month (slot 3);
  # pos4 = +3 days (slot 1, wrapped).
  expect_equal(out, as.POSIXct(c("2024-01-01", "2024-01-04", "2024-02-04", "2024-02-07"),
                               tz = attr(out, "tzone")))
})

test_that("idx_to_date errors for a multi_step_pattern calendar with posixct = FALSE", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), multi_step = msp, posixct = FALSE)
  expect_error(idx_to_date(cal, 1:3), "posixct = TRUE")
})

test_that("idx_calendar's existing numeric pattern mechanism is unaffected by multi_step_pattern support", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days",
                      pattern = c(1, 1, 1, 1, 3), pattern_start = 1L,
                      posixct = TRUE)
  expect_null(cal$multi_step)
  out <- idx_to_date(cal, 1:6)
  expect_equal(out, as.Date(c("2024-01-01", "2024-01-02", "2024-01-03",
                              "2024-01-04", "2024-01-05", "2024-01-08")))
})


test_that("xts_to_idx round-trips a real dataset's anchor date through idx_to_date/idx_to_pos", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  expect_true(is_idx_series(conv$series))
  expect_true(is_idx_calendar(conv$calendar))
  
  first_date <- zoo::index(gauteng)[1]
  expect_equal(idx_to_date(conv$calendar, conv$series$start), first_date)
  expect_equal(idx_to_pos(conv$calendar, first_date), conv$series$start)
})

test_that("idx_axis_opts constructs with defaults and validates mode", {
  ax <- idx_axis_opts()
  expect_s3_class(ax, "idx_axis_opts")
  expect_equal(ax$mode, "auto")
  expect_false(ax$info_box)
  expect_null(ax$pattern_n)
  
  ax2 <- idx_axis_opts(mode = "steps")
  expect_equal(ax2$mode, "steps")
  
  expect_error(idx_axis_opts(mode = "not_a_mode"))
})

test_that("idx_axis_opts validates info_box and pattern_n", {
  expect_error(idx_axis_opts(info_box = "TRUE"), "info_box")
  expect_error(idx_axis_opts(info_box = NA), "info_box")
  expect_error(idx_axis_opts(info_box = c(TRUE, FALSE)), "info_box")
  
  expect_error(idx_axis_opts(pattern_n = 0), "pattern_n")
  expect_error(idx_axis_opts(pattern_n = -1), "pattern_n")
  expect_error(idx_axis_opts(pattern_n = c(1, 2)), "pattern_n")
  expect_error(idx_axis_opts(pattern_n = "a"), "pattern_n")
  
  ax <- idx_axis_opts(mode = "time_since", info_box = TRUE, pattern_n = 3)
  expect_true(ax$info_box)
  expect_equal(ax$pattern_n, 3L)
  expect_type(ax$pattern_n, "integer")
})

test_that("idx_axis_opts modes 'steps', 'time_since' and 'date' require a calendar downstream", {
  set.seed(1)
  Y <- idx_series(cumsum(rpois(60, 8)) + 1, start = 1L)
  model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 50)
  res <- estimate(model)
  
  expect_error(plot_gy_ci(res, axis = idx_axis_opts(mode = "date")), "requires a calendar")
  expect_error(plot_gy_ci(res, axis = idx_axis_opts(mode = "steps")), "requires a calendar")
  expect_error(plot_gy_ci(res, axis = idx_axis_opts(mode = "time_since")), "requires a calendar")
  expect_no_error(plot_gy_ci(res, axis = idx_axis_opts(mode = "position")))
})
## Non-posixct mode: calendar core (idx_to_date, idx_to_pos, idx_offset_to_pos, idx_calendar_offset, idx_resolve_axis) ----

test_that("a non-posixct calendar can be built with a Date anchor and used end-to-end", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  
  # Same anchor date as xts_to_idx would use, but posixct = FALSE.
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  
  expect_false(cal_np$posixct)
  # idx_to_date falls back to plain Date + offset arithmetic.
  expect_equal(idx_to_date(cal_np, 1:5),
               zoo::index(gauteng)[1] + 0:4)
})

test_that("a non-posixct calendar can be built with a plain numeric anchor", {
  cal_np <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                         unit = "picoseconds", posixct = FALSE)
  expect_false(cal_np$posixct)
  expect_equal(idx_to_date(cal_np, 1:4), c(0, 2.5, 5, 7.5))
})

test_that("idx_to_pos rejects a non-posixct calendar even with a Date anchor", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  expect_error(idx_to_pos(cal_np, "2024-01-10"), "not a calendar-anchored")
})

test_that("idx_offset_to_pos is the correct inverse lookup for non-posixct calendars", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  # idx_to_date's plain-arithmetic branch fires regardless of anchor class
  # when posixct = FALSE, so idx_offset_to_pos is the correct inverse.
  d <- idx_to_date(cal_np, 1:10)
  expect_equal(idx_offset_to_pos(cal_np, d[7]), 7L)
  
  cal_ps <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                         unit = "picoseconds", posixct = FALSE)
  expect_equal(idx_offset_to_pos(cal_ps, 12.5), 6L)
})

test_that("non-posixct calendars respect a repeating pattern via idx_calendar_offset", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days",
                         pattern = c(1, 1, 1, 1, 3), pattern_start = 1L,
                         posixct = FALSE)
  expect_equal(idx_calendar_offset(cal_np, 5), 4)
  expect_equal(idx_calendar_offset(cal_np, 6), 7)
  out <- idx_to_date(cal_np, 1:6)
  expect_equal(out, as.Date("2024-01-01") + c(0, 1, 2, 3, 4, 7))
})

test_that("non-posixct calendars fall back to plain arithmetic for month/quarter/year units (no calendar-aware stepping)", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-31"), anchor_pos = 1L,
                         amount = 1, unit = "months", posixct = FALSE)
  out <- idx_to_date(cal_np, 1:3)
  # Plain arithmetic (anchor + offset * amount), not calendar-aware
  # seq.Date(by = "month") stepping.
  expect_equal(out, as.Date("2024-01-31") + 0:2)
  expect_false(isTRUE(all.equal(out, seq(as.Date("2024-01-31"), by = "month", length.out = 3))))
})

test_that("idx_to_date: posixct=TRUE but a numeric (non-Date/POSIXct) anchor still falls back to plain arithmetic", {
  # posixct = TRUE alone is not sufficient without a real calendar anchor.
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 1,
                      unit = "days", posixct = TRUE)
  expect_equal(idx_to_date(cal, 1:3), c(0, 1, 2))
})

test_that("idx_to_date: posixct=TRUE + Date anchor + non-standard unit still falls back to plain arithmetic", {
  # A non-standard unit has no recognised by_unit, so calendar-aware
  # stepping is skipped even with posixct = TRUE and a Date anchor.
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 10, unit = "fortnights", posixct = TRUE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, as.Date("2024-01-01") + c(0, 10, 20))
})

test_that("idx_to_date: posixct=TRUE + POSIXct anchor + standard sub-day unit uses calendar-aware seq.POSIXt stepping (contrast case)", {
  anchor <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  cal <- idx_calendar(anchor = anchor, anchor_pos = 1L,
                      amount = 1, unit = "hours", posixct = TRUE)
  out <- idx_to_date(cal, 1:5)
  expected <- seq(anchor, by = "hour", length.out = 5)
  # Compare by underlying instant, not class/attributes: idx_to_date()
  # does not preserve POSIXct class or tzone through coercion.
  expect_equal(as.numeric(out), as.numeric(expected))
})

test_that("idx_to_date: calendar-aware stepping covers the five Date-compatible standard units for a Date anchor (contrast case)", {
  units <- c("days", "weeks", "months", "quarters", "years")
  for (u in units) {
    cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                        amount = 1, unit = u, posixct = TRUE)
    out <- idx_to_date(cal, 1:3)
    expected <- seq(as.Date("2024-01-01"), by = sub("s$", "", u), length.out = 3)
    expect_equal(out, expected, info = paste("unit =", u))
  }
})

test_that("idx_to_date: calendar-aware stepping covers the sub-day standard units for a POSIXct anchor (contrast case)", {
  anchor <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  # seq.POSIXt wants abbreviated "sec"/"min" forms, matching
  # idx_calendar_by_unit()'s mapping.
  by_strings <- c(seconds = "sec", minutes = "min", hours = "hour")
  for (u in names(by_strings)) {
    cal <- idx_calendar(anchor = anchor, anchor_pos = 1L,
                        amount = 1, unit = u, posixct = TRUE)
    out <- idx_to_date(cal, 1:3)
    expected <- seq(anchor, by = by_strings[[u]], length.out = 3)
    # Compare by underlying instant, not strict object equality.
    expect_equal(as.numeric(out), as.numeric(expected), info = paste("unit =", u))
  }
})

test_that("idx_to_date: unit matching is case-insensitive and singular/plural-tolerant (calendar-aware branch)", {
  cal_plural <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                             amount = 1, unit = "Days", posixct = TRUE)
  cal_singular <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                               amount = 1, unit = "DAY", posixct = TRUE)
  expected <- seq(as.Date("2024-01-01"), by = "day", length.out = 3)
  expect_equal(idx_to_date(cal_plural, 1:3), expected)
  expect_equal(idx_to_date(cal_singular, 1:3), expected)
})

test_that("idx_to_date: posixct=FALSE + POSIXct anchor + standard unit -> plain arithmetic, preserves POSIXct class", {
  anchor <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  cal <- idx_calendar(anchor = anchor, anchor_pos = 1L,
                      amount = 3600, unit = "seconds", posixct = FALSE)
  out <- idx_to_date(cal, 1:3)
  expect_equal(out, anchor + c(0, 3600, 7200))
  expect_s3_class(out, "POSIXct")
})

test_that("idx_to_date: plain-arithmetic fallback is vectorized and handles negative/out-of-order positions under posixct = FALSE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 10L,
                      amount = 1, unit = "days", posixct = FALSE)
  out <- idx_to_date(cal, c(15, 5, 10, -5))
  expect_equal(out, as.Date("2024-01-01") + c(5, -5, 0, -15))
})

test_that("idx_to_date: non-integer step check is bypassed entirely when posixct = FALSE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 2.5, unit = "days", posixct = FALSE)
  expect_no_error(out <- idx_to_date(cal, 2))
  expect_equal(out, as.Date("2024-01-01") + 2.5)
})

test_that("idx_to_pos errors when posixct = FALSE even though the anchor is a valid POSIXct", {
  anchor <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  cal <- idx_calendar(anchor = anchor, anchor_pos = 1L,
                      amount = 1, unit = "hours", posixct = FALSE)
  expect_error(idx_to_pos(cal, "2024-01-01 05:00:00"), "not a calendar-anchored")
})

test_that("idx_to_pos error message points the user at idx_offset_to_pos() for non-calendar-anchored calendars", {
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 1, unit = "days")
  expect_error(idx_to_pos(cal, "2024-01-10"), "idx_offset_to_pos")
})

test_that("idx_offset_to_pos works identically whether posixct is TRUE or FALSE (posixct-agnostic by design)", {
  cal_false <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                            unit = "picoseconds", posixct = FALSE)
  cal_true <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
                           unit = "picoseconds", posixct = TRUE)
  expect_equal(idx_offset_to_pos(cal_false, 12.5), 6L)
  expect_equal(idx_offset_to_pos(cal_true, 12.5), 6L)
})

test_that("idx_offset_to_pos works with a Date anchor even though idx_to_pos() would reject a non-posixct one", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = FALSE)
  expect_equal(idx_offset_to_pos(cal, as.Date("2024-01-05")), 5L)
})

test_that("idx_resolve_axis: auto resolves to position when no calendar is supplied", {
  ax <- idx_resolve_axis(NULL, NULL)
  expect_equal(ax$mode, "position")
})

test_that("idx_resolve_axis: 'no calendar' and 'posixct=FALSE calendar' give the identical auto resolution", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = FALSE)
  expect_equal(idx_resolve_axis(NULL, NULL)$mode, idx_resolve_axis(NULL, cal)$mode)
})

test_that("idx_resolve_axis: explicit mode='position' works with a posixct=FALSE calendar, a posixct=TRUE calendar, or no calendar", {
  cal_false <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  cal_true  <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = TRUE)
  ax_opts <- idx_axis_opts(mode = "position")
  
  expect_equal(idx_resolve_axis(ax_opts, NULL)$mode, "position")
  expect_equal(idx_resolve_axis(ax_opts, cal_false)$mode, "position")
  expect_equal(idx_resolve_axis(ax_opts, cal_true)$mode, "position")
})

test_that("idx_resolve_axis: explicit mode='date' succeeds even with a numeric (non-Date/POSIXct) posixct=FALSE anchor", {
  # idx_resolve_axis() only checks is_idx_calendar(); the fallback to
  # plain arithmetic happens downstream in idx_to_date().
  cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 1,
                      unit = "widgets", posixct = FALSE)
  ax <- idx_resolve_axis(idx_axis_opts(mode = "date"), cal)
  expect_equal(ax$mode, "date")
})

test_that("'auto' axis mode resolves to 'position' (not 'date') under a non-posixct calendar", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  resolved <- idx_resolve_axis(idx_axis_opts(mode = "auto"), cal_np)
  expect_equal(resolved$mode, "position")
})

test_that("'date' axis mode is still permitted (calendar is present) but plots a plain-arithmetic x value under a non-posixct calendar", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  cal_np <- idx_calendar(anchor = zoo::index(gauteng)[1], anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  end.pos <- idx_to_pos(conv$calendar, as.Date("2020-07-20"))
  model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005, end = end.pos, calendar = cal_np)
  res <- estimate(model)
  
  # "date" mode only requires a calendar object, not posixct = TRUE.
  expect_no_error(plot_forecast(res, axis = idx_axis_opts(mode = "date")))
  
  df <- idx_series_df(conv$series, cal_np, idx_axis_opts(mode = "date"))
  expect_equal(df$x, idx_to_date(cal_np, idx_positions(conv$series)))
})

test_that("xts_to_idx() always sets posixct = TRUE (contrast case documenting why manual non-posixct calendars are needed to exercise this mode)", {
  data(gauteng, package = "tsgc")
  conv <- xts_to_idx(gauteng)
  expect_true(conv$calendar$posixct)
})
## Non-posixct mode: idx_series_df / idx_x_lab / idx_x_scale / idx_info_box_caption / idx_add_info_box ----

test_that("idx_series_df: mode='position' returns raw positions regardless of calendar/posixct", {
  x <- idx_series(c(10, 20, 30), start = 5L)
  df_no_cal <- idx_series_df(x, NULL, idx_axis_opts(mode = "position"))
  expect_equal(df_no_cal$x, 5:7)
  
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 5L,
                      amount = 1, unit = "days", posixct = FALSE)
  df_with_cal <- idx_series_df(x, cal, idx_axis_opts(mode = "position"))
  expect_equal(df_with_cal$x, 5:7)
})

test_that("idx_series_df: mode='steps' computes position - anchor_pos exactly, under posixct = FALSE", {
  x <- idx_series(c(10, 20, 30), start = 5L)
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 3L,
                      amount = 1, unit = "days", posixct = FALSE)
  df <- idx_series_df(x, cal, idx_axis_opts(mode = "steps"))
  expect_equal(df$x, c(2, 3, 4))
})

test_that("idx_series_df: mode='time_since' applies the pattern-weighted offset x amount exactly, under posixct = FALSE", {
  x <- idx_series(c(10, 20, 30, 40, 50), start = 1L)
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 2, unit = "days",
                      pattern = c(1, 1, 1, 1, 3), pattern_start = 1L,
                      posixct = FALSE)
  df <- idx_series_df(x, cal, idx_axis_opts(mode = "time_since"))
  # idx_calendar_offset at positions 1:5 (Mon-Fri, all within first cycle):
  # offsets 0, 1, 2, 3, 4, times amount = 2.
  expect_equal(df$x, c(0, 1, 2, 3, 4) * 2)
})

test_that("idx_series_df: mode='date' under posixct=FALSE returns the plain-arithmetic date exactly, still Date-classed", {
  x <- idx_series(c(10, 20, 30), start = 1L)
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = FALSE)
  df <- idx_series_df(x, cal, idx_axis_opts(mode = "date"))
  expect_equal(df$x, as.Date("2024-01-01") + 0:2)
  expect_s3_class(df$x, "Date")
})


test_that("idx_series_df: mode='auto' falls back to plain positions for a posixct=FALSE calendar", {
  x <- idx_series(c(10, 20, 30), start = 5L)
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                      amount = 1, unit = "days", posixct = FALSE)
  df <- idx_series_df(x, cal, NULL)
  expect_equal(df$x, 5:7)
})

test_that("idx_x_lab: 'position' label is always 'Position', with or without a calendar", {
  expect_equal(idx_x_lab(NULL, idx_axis_opts(mode = "position")), "Position")
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  expect_equal(idx_x_lab(cal, idx_axis_opts(mode = "position")), "Position")
})

test_that("idx_x_lab: 'steps' label falls back to generic 'anchor' wording when anchor_name is unset, under posixct = FALSE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  expect_equal(idx_x_lab(cal, idx_axis_opts(mode = "steps")), "Steps from anchor")
})

test_that("idx_x_lab: 'steps' label uses anchor_name when set, under posixct = FALSE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE,
                      anchor_name = "product launch")
  expect_equal(idx_x_lab(cal, idx_axis_opts(mode = "steps")), "Steps from product launch")
})

test_that("idx_x_lab: 'time_since' label includes calendar$unit, under posixct = FALSE", {
  cal <- idx_calendar(anchor = 0, unit = "picoseconds", posixct = FALSE,
                      anchor_name = "t-zero")
  expect_equal(idx_x_lab(cal, idx_axis_opts(mode = "time_since")), "Time since t-zero (picoseconds)")
})

test_that("idx_x_lab: 'date' label is always 'Date', irrespective of posixct", {
  cal_false <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  cal_true  <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = TRUE)
  expect_equal(idx_x_lab(cal_false, idx_axis_opts(mode = "date")), "Date")
  expect_equal(idx_x_lab(cal_true, idx_axis_opts(mode = "date")), "Date")
})

test_that("idx_x_scale: non-date modes always return a continuous scale, regardless of posixct", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  sc <- idx_x_scale(cal, idx_axis_opts(mode = "position"))
  expect_true(inherits(sc, "Scale"))
  expect_false(grepl("Date|Datetime", paste(class(sc), collapse = " ")))
})

test_that("idx_x_scale: explicit mode='date' on a posixct=FALSE calendar chooses scale_x_date/scale_x_datetime by anchor class, not by the posixct flag", {
  cal_date <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  cal_psx  <- idx_calendar(anchor = as.POSIXct("2024-01-01", tz = "UTC"), posixct = FALSE)
  
  sc_date <- idx_x_scale(cal_date, idx_axis_opts(mode = "date"))
  sc_psx  <- idx_x_scale(cal_psx, idx_axis_opts(mode = "date"))
  
  expect_true(grepl("Date", paste(class(sc_date), collapse = " ")))
  expect_true(grepl("Datetime", paste(class(sc_psx), collapse = " ")))
})

test_that("idx_x_scale: explicit mode='date' does not error for a numeric (non-Date/POSIXct) posixct=FALSE anchor", {
  cal <- idx_calendar(anchor = 0, posixct = FALSE)
  expect_no_error(idx_x_scale(cal, idx_axis_opts(mode = "date")))
})

test_that("idx_info_box_caption: formats a multi_step_pattern calendar by describing each step slot", {
  msp <- multi_step_pattern(idx_step(days = 3), idx_step(months = 1))
  cal <- idx_calendar_multi_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                                 multi_step = msp, posixct = TRUE,
                                 anchor_name = "launch")
  cap <- idx_info_box_caption(cal, idx_axis_opts(mode = "steps", info_box = TRUE))
  expect_true(grepl("Anchor: launch", cap, fixed = TRUE))
  expect_true(grepl("multi_step_pattern \\(2 slots\\)", cap))
  expect_true(grepl("\\[1\\] 3 days", cap))
  expect_true(grepl("\\[2\\] 1 months", cap))
})

test_that("idx_info_box_caption: formats an idx_calendar_step (compound step, NA amount) calendar via idx_step_label", {
  cal <- idx_calendar_step(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                           step = idx_step(months = 1, days = 5), posixct = TRUE)
  cap <- idx_info_box_caption(cal, idx_axis_opts(mode = "steps", info_box = TRUE))
  expect_true(grepl("Step: 1 months, 5 days", cap))
})

test_that("idx_info_box_caption: NULL when info_box is FALSE/unset or there is no calendar, regardless of posixct", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE)
  expect_null(idx_info_box_caption(cal, NULL))
  expect_null(idx_info_box_caption(cal, idx_axis_opts(mode = "position", info_box = FALSE)))
  expect_null(idx_info_box_caption(NULL, idx_axis_opts(mode = "position", info_box = TRUE)))
})

test_that("idx_info_box_caption: formats a numeric (non-Date) anchor correctly under posixct = FALSE", {
  cal <- idx_calendar(anchor = 0, amount = 2.5, unit = "picoseconds", posixct = FALSE)
  cap <- idx_info_box_caption(cal, idx_axis_opts(mode = "position", info_box = TRUE))
  expect_true(grepl("Anchor: 0", cap))
  expect_true(grepl("Step: 2.5 picoseconds", cap))
  expect_true(grepl("Pattern: 1", cap))
})

test_that("idx_info_box_caption: includes anchor_name when set, under posixct = FALSE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE,
                      anchor_name = "launch")
  cap <- idx_info_box_caption(cal, idx_axis_opts(mode = "position", info_box = TRUE))
  expect_true(grepl("Anchor: launch \\(2024-01-01\\)", cap))
})

test_that("idx_info_box_caption: pattern_n truncates a long pattern under posixct = FALSE", {
  cal <- idx_calendar(anchor = as.Date("2024-01-01"), posixct = FALSE,
                      pattern = c(1, 1, 1, 1, 1, 1, 3), pattern_start = 1L)
  cap_full <- idx_info_box_caption(cal, idx_axis_opts(mode = "position", info_box = TRUE))
  cap_trunc <- idx_info_box_caption(cal, idx_axis_opts(mode = "position", info_box = TRUE, pattern_n = 3))
  expect_true(grepl("Pattern: 1, 1, 1, 1, 1, 1, 3", cap_full))
  expect_true(grepl("Pattern: 1, 1, 1, \\.\\.\\.", cap_trunc))
})

test_that("idx_add_info_box appends to (not overwrites) an existing caption, under posixct = FALSE", {
  cal <- idx_calendar(anchor = 0, posixct = FALSE)
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_line() +
    ggplot2::labs(caption = "existing caption")
  p2 <- idx_add_info_box(p, cal, idx_axis_opts(mode = "position", info_box = TRUE))
  expect_true(grepl("existing caption", p2$labels$caption))
  expect_true(grepl("Anchor: 0", p2$labels$caption))
})

test_that("idx_x_lab() and idx_info_box_caption() behave sensibly for non-posixct calendars", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE,
                         anchor_name = "launch")
  expect_equal(idx_x_lab(cal_np, idx_axis_opts(mode = "steps")), "Steps from launch")
  expect_equal(idx_x_lab(cal_np, idx_axis_opts(mode = "time_since")), "Time since launch (days)")
  expect_equal(idx_x_lab(cal_np, idx_axis_opts(mode = "position")), "Position")
  
  cap <- idx_info_box_caption(cal_np, idx_axis_opts(mode = "steps", info_box = TRUE))
  expect_true(grepl("Anchor: launch", cap, fixed = TRUE))
})

test_that("idx_x_scale() falls back to a continuous scale (not scale_x_date/datetime) for non-'date' modes with a non-posixct calendar", {
  cal_np <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
                         amount = 1, unit = "days", posixct = FALSE)
  sc <- idx_x_scale(cal_np, idx_axis_opts(mode = "steps"))
  expect_true(inherits(sc, "ScaleContinuousPosition"))
})
test_that("idx_to_date preserves Date/POSIXct class for idx_calendar_multi_step calendars", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-01-01"),
    multi_step = multi_step_pattern(
      idx_step(days = 3), idx_step(days = 3), idx_step(months = 1)
    ),
    posixct = TRUE
  )
  out <- idx_to_date(cal_multi, 1:6)
  expect_true(inherits(out, "POSIXct"))
  expect_false(is.numeric(out) && !inherits(out, c("Date", "POSIXct")))
})

test_that("idx_to_date returns exact expected dates for idx_calendar_multi_step", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-01-01"),
    multi_step = multi_step_pattern(
      idx_step(days = 3), idx_step(days = 3), idx_step(months = 1)
    ),
    posixct = TRUE
  )
  out <- idx_to_date(cal_multi, 1:6)
  expected <- as.POSIXct(c(
    "2024-01-01", "2024-01-04", "2024-01-07",
    "2024-02-07", "2024-02-10", "2024-02-13"
  ), tz = attr(out, "tzone"))
  expect_equal(as.numeric(out), as.numeric(expected))
})

test_that("idx_to_pos(idx_to_date(multi_step_cal, p)) == p (round trip)", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-01-01"),
    multi_step = multi_step_pattern(
      idx_step(days = 3), idx_step(days = 3), idx_step(months = 1)
    ),
    posixct = TRUE
  )
  p <- 1:12
  dts <- idx_to_date(cal_multi, p)
  back <- vapply(seq_along(p), function(i) idx_to_pos(cal_multi, dts[i]), integer(1))
  expect_equal(back, p)
})

test_that("idx_calendar_multi_step handles positions before and after the anchor", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-01-01"),
    anchor_pos = 5L,
    multi_step = multi_step_pattern(
      idx_step(days = 3), idx_step(days = 3), idx_step(months = 1)
    ),
    posixct = TRUE
  )
  out <- idx_to_date(cal_multi, 3:7)
  expect_true(inherits(out, "POSIXct"))
  expect_equal(idx_to_pos(cal_multi, out[1]), 3L)
  expect_equal(idx_to_pos(cal_multi, out[5]), 7L)
})

test_that("idx_calendar_multi_step handles month-end boundaries without losing class", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-01-30"),
    multi_step = multi_step_pattern(
      idx_step(months = 1), idx_step(days = 1)
    ),
    posixct = TRUE
  )
  out <- idx_to_date(cal_multi, 1:4)
  expect_true(inherits(out, "POSIXct"))
  expect_false(any(is.na(out)))
})

test_that("idx_calendar_multi_step handles leap-year boundaries without losing class", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-02-28"),
    multi_step = multi_step_pattern(
      idx_step(days = 1), idx_step(days = 1)
    ),
    posixct = TRUE
  )
  out <- idx_to_date(cal_multi, 1:4)
  expect_true(inherits(out, "POSIXct"))
  expect_equal(as.Date(out[3]), as.Date("2024-03-01"))
})

test_that("idx_calendar_multi_step handles a mixed-step (fixed + calendar-relative) cycle without losing class", {
  cal_multi <- idx_calendar_multi_step(
    anchor = as.Date("2024-01-01"),
    multi_step = multi_step_pattern(
      idx_step(weeks = 1), idx_step(quarters = 1), idx_step(days = 5)
    ),
    posixct = TRUE
  )
  out <- idx_to_date(cal_multi, 1:9)
  expect_true(inherits(out, "POSIXct"))
  expect_false(is.numeric(out))
  expect_true(all(diff(as.numeric(out)) > 0))
})

test_that("idx_to_pos exactly inverts numeric-pattern calendars", {
  bdays <- seq(as.Date("2024-01-01"), by = "day", length.out = 20)
  bdays <- bdays[!weekdays(bdays) %in% c("Saturday", "Sunday")]
  x <- xts::xts(seq_along(bdays), order.by = bdays)
  cal <- xts_to_idx(x)$calendar
  
  back <- vapply(seq_along(bdays), function(i) idx_to_pos(cal, bdays[i]), integer(1))
  expect_equal(back, seq_along(bdays))
})

test_that("idx_to_pos rejects dates between whole numeric-pattern steps", {
  cal <- idx_calendar(as.Date("2024-01-01"), amount = 2, unit = "days",
                      posixct = TRUE)
  expect_error(idx_to_pos(cal, "2024-01-02"), "whole step boundary")
})