utils::globalVariables(c("Date", "Rt", "lower", "upper", "forecast", "model", "x", "Centered.MA"))

#' @importFrom utils head tail
#' @importFrom stats setNames
NULL

#' @title Detect the shortest repeating step pattern in a calendar index
#'
#' @description Internal helper for \code{\link{xts_to_idx}}. Given a
#' sorted, unique vector of \code{Date}/\code{POSIXct}/\code{yearqtr}/
#' \code{yearmon} values, finds a single \code{amount}/\code{unit} pair -
#' searched across the full unit ladder \code{years}, \code{quarters},
#' \code{months}, \code{weeks}, \code{days}, \code{hours}, \code{minutes},
#' \code{seconds} - such that every gap between consecutive values is a
#' whole-number multiple of it, then finds the shortest repeating cycle of
#' those multiples (e.g. business days -> unit \code{"days"}, cycle
#' \code{c(1, 1, 1, 1, 3)}). \code{years}/\code{quarters}/\code{months} are
#' calendar-relative (variable length in days/seconds), so - for
#' \code{Date}/\code{POSIXct} indices - these are checked via calendar-aware
#' whole-month counting (mirroring \code{\link{idx_to_pos}}), not via a
#' fixed-duration ratio; \code{weeks} through \code{seconds} are
#' fixed-duration and checked via a plain ratio. The largest unit that
#' fits is preferred (e.g. a quarterly index is reported as
#' \code{amount = 1, unit = "quarters"}, not \code{amount = 91,
#' unit = "days"}), since that is what a person building the equivalent
#' \code{idx_calendar} by hand would write.
#'
#' This is deliberately conservative: it only recognises a single
#' \code{amount}/\code{unit} pair (scaled by an integer vector
#' \code{pattern}), not the more general compound \code{idx_step}/
#' \code{multi_step_pattern} constructors - those still require manual
#' construction via \code{\link{idx_calendar_step}}/
#' \code{\link{idx_calendar_multi_step}}. Returns \code{NULL} (rather than
#' erroring) if no such single-unit pattern fits.
#'
#' @param idx A sorted, unique vector of length >= 2, of class \code{Date},
#' \code{POSIXct}, \code{yearqtr}, or \code{yearmon}.
#' @param max_pattern_len Maximum repeating-cycle length to search for.
#' Defaults to \code{7} (enough for business-day-style weekly cycles);
#' raising this is rarely useful and slows detection on long, irregular
#' indices, since every candidate length up to it is checked.
#'
#' @returns \code{NULL}, or a list with elements \code{amount}, \code{unit},
#' \code{pattern} (integer vector) suitable for passing to
#' \code{\link{idx_calendar}}.
#' @keywords internal
#' @noRd
idx_detect_calendar_pattern <- function(idx, max_pattern_len = 7L) {
  if (length(idx) < 2) return(NULL)
  
  is_yearqtr <- inherits(idx, "yearqtr")
  is_yearmon <- inherits(idx, "yearmon")
  
  # Shortest repeating cycle in a vector of positive whole-unit gap multiples.
  find_pattern <- function(mult) {
    if (any(mult <= 0) || !isTRUE(all.equal(mult, round(mult)))) return(NULL)
    mult <- round(mult)
    m <- length(mult)
    cap <- min(max_pattern_len, m)
    for (plen in seq_len(cap)) {
      candidate <- mult[seq_len(plen)]
      full_cycles <- mult[seq_len(m - m %% plen)]
      reps <- matrix(full_cycles, nrow = plen)
      ok <- all(apply(reps, 2, function(col) isTRUE(all.equal(col, candidate))))
      if (ok) {
        remainder <- m %% plen
        if (remainder > 0) {
          tail_part <- mult[(m - remainder + 1L):m]
          ok <- isTRUE(all.equal(tail_part, candidate[seq_len(remainder)]))
        }
        if (ok) return(as.integer(candidate))
      }
    }
    NULL
  }
  
  if (is_yearqtr || is_yearmon) {
    # yearqtr/yearmon are fractional-year numerics; try the native unit
    # (quarters/months), then the coarser years unit.
    base_step <- if (is_yearqtr) 0.25 else (1 / 12)
    native_unit <- if (is_yearqtr) "quarters" else "months"
    gaps <- diff(as.numeric(idx))
    
    years_mult <- gaps / 1
    pat <- find_pattern(years_mult)
    if (!is.null(pat)) return(list(amount = 1, unit = "years", pattern = pat))
    
    native_mult <- gaps / base_step
    pat <- find_pattern(native_mult)
    if (!is.null(pat)) return(list(amount = 1, unit = native_unit, pattern = pat))
    
    return(NULL)
  }
  
  # Date/POSIXct: try calendar-relative units first (largest to smallest).
  # Only attempted when idx has no time-of-day component to lose.
  no_time_of_day <- if (inherits(idx, "POSIXct")) {
    all(format(idx, "%H:%M:%S") == "00:00:00")
  } else {
    TRUE
  }
  if (no_time_of_day) {
    ymd <- as.Date(idx)
    ym <- as.integer(format(ymd, "%Y")) * 12L + (as.integer(format(ymd, "%m")) - 1L)
    months_gaps <- diff(ym)
    # Only calendar-relative if every value falls on the first of its month.
    on_month_start <- all(as.integer(format(ymd, "%d")) == 1L)
    if (on_month_start && all(months_gaps > 0)) {
      for (per_year in c(1, 4, 12)) {
        unit <- switch(as.character(per_year), "1" = "years", "4" = "quarters", "12" = "months")
        step_months <- 12 / per_year
        pat <- find_pattern(months_gaps / step_months)
        if (!is.null(pat)) return(list(amount = 1, unit = unit, pattern = pat))
      }
    }
  }
  
  # Fixed-duration units: work in whole seconds so Date and POSIXct share
  # one code path (a Date's "seconds" are just days * 86400).
  secs <- if (inherits(idx, "POSIXct")) as.numeric(idx) else as.numeric(idx) * 86400
  gap_secs <- diff(secs)
  if (any(gap_secs <= 0)) return(NULL)
  
  candidate_units <- list(weeks = 86400 * 7, days = 86400, hours = 3600,
                          minutes = 60, seconds = 1)
  for (nm in names(candidate_units)) {
    u <- candidate_units[[nm]]
    ratios <- gap_secs / u
    if (all(abs(ratios - round(ratios)) < 1e-6)) {
      pat <- find_pattern(ratios)
      if (!is.null(pat)) return(list(amount = 1, unit = nm, pattern = pat))
    }
  }
  NULL
}

#' @title Convert an \code{xts}/\code{zoo} object to an \code{idx_series}
#'
#' @description Converts a calendar-indexed \code{xts} or \code{zoo} object
#' into an \code{\link{idx_series}}, together with an
#' \code{\link{idx_calendar}} describing how to translate back to calendar
#' time. This is the standard entry point for bringing calendar-indexed
#' data (the user's own \code{xts}/\code{zoo}/data frame data, or one of
#' the package's bundled datasets) into the \code{idx_series}-based
#' analysis functions.
#'
#' By default (\code{detect = TRUE}), the step size and any repeating gap
#' pattern are auto-detected from \code{x}'s index via
#' \code{idx_detect_calendar_pattern} (searching the full unit
#' ladder years/quarters/months/weeks/days/hours/minutes/seconds), so
#' callers do not need to work out \code{amount}/\code{unit}/\code{pattern}
#' themselves for common cases: a plain daily/weekly/etc. index, a
#' business-day index (weekends skipped), a monthly/quarterly/yearly
#' \code{Date} index, or a \code{yearqtr}/\code{yearmon} index. This only
#' recognises a single \code{amount}/\code{unit} pair (optionally repeated
#' via an integer vector \code{pattern}); it does not attempt to detect
#' compound (\code{idx_step}) or heterogeneous (\code{multi_step_pattern})
#' calendars - build those manually via \code{\link{idx_calendar_step}}/
#' \code{\link{idx_calendar_multi_step}} if needed. If detection fails
#' (the index has no such regular pattern, e.g. genuinely irregular gaps),
#' \code{xts_to_idx} raises an error explaining why, rather than silently
#' falling back to an incorrect calendar; pass \code{detect = FALSE} (and
#' set \code{amount}/\code{unit} yourself) for indices detection can't
#' handle.
#'
#' @param x An \code{xts} or \code{zoo} object with a \code{Date},
#' \code{POSIXct}, \code{yearqtr}, or \code{yearmon} index of at least 2
#' rows (needed to detect a step size), and no duplicate index values.
#' @param start.pos Integer position that the first row of \code{x} should
#' be assigned. Use this to align \code{x} with another, already-converted
#' \code{idx_series} that starts at a different calendar date - e.g.
#' \code{start.pos = idx_to_pos(other_cal, zoo::index(x)[1])}. Defaults to
#' \code{1L}.
#' @param detect A single logical. If \code{TRUE} (default), auto-detect
#' \code{amount}/\code{unit}/\code{pattern} from \code{x}'s index (see
#' Description). If \code{FALSE}, use \code{amount}/\code{unit} as
#' supplied (previous behaviour; defaults to \code{amount = 1,
#' unit = "days"} to match this function's original hardcoded assumption).
#' @param amount,unit Used only when \code{detect = FALSE}; passed straight
#' through to \code{\link{idx_calendar}}.
#'
#' @returns A list with two elements: \code{series}, an \code{idx_series}
#' holding \code{x}'s values, and \code{calendar}, an \code{idx_calendar}
#' anchoring \code{series}'s positions to \code{x}'s original dates.
#'
#' @examples
#' # Daily data - detected automatically.
#' x <- xts::xts(cumsum(rpois(30, 5)) + 1, order.by = Sys.Date() - 29:0)
#' conv <- xts_to_idx(x)
#' conv$series
#' conv$calendar
#'
#' # Business-day data (weekends skipped) - pattern detected automatically.
#' bdays <- seq(as.Date("2024-01-01"), by = "day", length.out = 40)
#' bdays <- bdays[!weekdays(bdays) %in% c("Saturday", "Sunday")]
#' xb <- xts::xts(cumsum(rpois(length(bdays), 5)) + 1, order.by = bdays)
#' conv_b <- xts_to_idx(xb)
#' conv_b$calendar
#'
#' # Monthly Date index - detected as amount = 1, unit = "months" (not days).
#' months <- seq(as.Date("2024-01-01"), by = "month", length.out = 24)
#' xm <- xts::xts(cumsum(rpois(24, 5)) + 1, order.by = months)
#' conv_m <- xts_to_idx(xm)
#' conv_m$calendar
#'
#' @importFrom xts xts
#' @importFrom zoo index coredata
#'
#' @export
xts_to_idx <- function(x, start.pos = 1L, detect = TRUE, amount = 1, unit = "days") {
  idx <- index(x)
  if (anyDuplicated(idx)) {
    stop("xts_to_idx: x's index contains duplicate values.")
  }
  ord <- order(idx)
  if (is.unsorted(ord)) {
    idx <- idx[ord]
    x <- x[ord]
  }
  
  if (isTRUE(detect)) {
    if (length(idx) < 2) {
      stop("xts_to_idx: detect = TRUE requires at least 2 rows to detect a ",
           "step size; pass detect = FALSE with an explicit amount/unit ",
           "for single-row data.")
    }
    if (!inherits(idx, c("Date", "POSIXct", "yearqtr", "yearmon"))) {
      stop("xts_to_idx: detect = TRUE requires a Date/POSIXct/yearqtr/",
           "yearmon index; got class ", paste(class(idx), collapse = "/"), ".")
    }
    detected <- idx_detect_calendar_pattern(idx)
    if (is.null(detected)) {
      stop("xts_to_idx: could not detect a regular single-unit step ",
           "pattern in x's index (gaps are irregular, or the repeating ",
           "cycle is longer than idx_detect_calendar_pattern's default ",
           "search length). Pass detect = FALSE and build the ",
           "idx_calendar yourself (idx_calendar_step()/",
           "idx_calendar_multi_step() for compound/heterogeneous steps).")
    }
    return(list(
      series = idx_series(coredata(x), start = start.pos),
      calendar = idx_calendar(
        anchor = idx[1],
        anchor_pos = start.pos,
        amount = detected$amount, unit = detected$unit,
        pattern = detected$pattern,
        posixct = TRUE
      )
    ))
  }
  
  list(
    series = idx_series(coredata(x), start = start.pos),
    calendar = idx_calendar(
      anchor = idx[1],
      anchor_pos = start.pos,
      amount = amount, unit = unit,
      posixct = TRUE
    )
  )
}

#' @title Translate a calendar date to an \code{idx_series} integer position
#'
#' @description Inverse of \code{\link{idx_to_date}}: converts a calendar
#' date into the integer position (relative to \code{cal}) that corresponds
#' to it. Intended for locating estimation start/end dates, or for aligning
#' a second series to an already-converted \code{idx_series} (see
#' \code{\link{xts_to_idx}}). Requires a calendar-anchored \code{cal}
#' (\code{cal$posixct = TRUE} and a \code{Date}/\code{POSIXct}/
#' \code{yearqtr}/\code{yearmon} \code{cal$anchor}); use
#' \code{\link{idx_offset_to_pos}} instead for non-calendar (e.g.
#' sub-second/arbitrary numeric) calendars.
#'
#' \code{cal} may be built via any of the three \code{idx_calendar}
#' constructors:
#' \itemize{
#'   \item The primary \code{\link{idx_calendar}} constructor (a single
#'   \code{amount}/\code{unit} pair). For the calendar-relative units
#'   \code{"months"}, \code{"quarters"}, \code{"years"}, the position is
#'   found by counting whole calendar steps between \code{cal$anchor} and
#'   \code{date} (the inverse of the calendar-aware stepping used by
#'   \code{\link{idx_to_date}}), rather than via a fixed-length difference -
#'   a quarter is not a fixed number of days. For the fixed-duration units
#'   \code{"seconds"}, \code{"minutes"}, \code{"hours"}, \code{"days"},
#'   \code{"weeks"}, positions are computed as a plain time difference
#'   scaled by \code{cal$amount}; this is exact for these units. A repeating
#'   numeric \code{pattern} is respected. Because patterned offsets are
#'   piecewise rather than a single amount that can be divided in closed
#'   form, the inverse is found by an exact binary search over integer
#'   positions. Dates that do not land on a whole step boundary are rejected
#'   rather than rounded.
#'   \item \code{\link{idx_calendar_step}} (a compound \code{idx_step}) and
#'   \code{\link{idx_calendar_multi_step}} (a \code{multi_step_pattern})
#'   calendars have no single \code{amount}/\code{unit} pair to invert in
#'   closed form in general, so \code{date} is instead located by search:
#'   \code{idx_calendar_step} calendars repeat one compound \code{idx_step}
#'   every position, and \code{\link{idx_step_add}} is strictly monotonic
#'   in the step count, so this uses an O(log distance) binary search;
#'   \code{idx_calendar_multi_step} calendars cycle through heterogeneous
#'   steps that cannot be searched this way, so this walks one step at a
#'   time from \code{cal$anchor_pos} instead, an O(distance) linear walk.
#'   Both are exact, and both ignore any numeric \code{pattern} (a
#'   \code{multi_step_pattern} calendar has no numeric \code{pattern} to
#'   ignore in the first place - it uses \code{pattern_start} only as the
#'   \code{multi_step} slot index).
#' }
#'
#' A \code{yearqtr}/\code{yearmon} \code{cal$anchor} is only valid together
#' with \code{cal$unit} of \code{"quarters"}/\code{"months"} respectively
#' (matching \code{\link{idx_step_add}}'s restriction that these anchors
#' have no sub-quarter/sub-month resolution); \code{date} is coerced to the
#' same class as \code{cal$anchor} and whole-step boundaries are found via
#' fractional-year arithmetic, consistent with \code{\link{idx_to_date}}.
#'
#' @param cal An \code{idx_calendar} object with \code{posixct = TRUE} and a
#' \code{Date}/\code{POSIXct}/\code{yearqtr}/\code{yearmon} \code{anchor}.
#' @param date A single date/time. For \code{cal$unit} of \code{"months"},
#' \code{"quarters"}, or \code{"years"}, coerced via \code{as.Date} (or to
#' \code{yearqtr}/\code{yearmon}, matching \code{cal$anchor}). For all
#' other units, coerced via \code{as.POSIXct} if \code{cal$anchor} is
#' \code{POSIXct}, or \code{as.Date} otherwise - so a character string with
#' a time-of-day component (e.g. \code{"2024-01-01 13:30:00"}) is preserved
#' where relevant (sub-day units).
#'
#' @returns A single integer: the \code{idx_series} position corresponding
#' to \code{date}.
#'
#' @examples
#' cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
#'                      amount = 1, unit = "days", posixct = TRUE)
#' idx_to_pos(cal, "2024-01-10")
#'
#' # A multi_step_pattern calendar (walked, not closed-form).
#' cal_ms <- idx_calendar_multi_step(
#'   anchor = as.Date("2024-01-01"),
#'   multi_step = multi_step_pattern(idx_step(days = 3), idx_step(days = 3),
#'                                    idx_step(months = 1)),
#'   posixct = TRUE
#' )
#' idx_to_pos(cal_ms, "2024-02-07")
#'
#' @export
idx_to_pos <- function(cal, date) {
  is_yearqtr_anchor <- inherits(cal$anchor, "yearqtr")
  is_yearmon_anchor <- inherits(cal$anchor, "yearmon")
  is_date_posix_anchor <- inherits(cal$anchor, "Date") || inherits(cal$anchor, "POSIXct")
  is_calendar_anchor <- is_date_posix_anchor || is_yearqtr_anchor || is_yearmon_anchor
  
  if (!isTRUE(cal$posixct) || !is_calendar_anchor) {
    stop("idx_to_pos: cal is not a calendar-anchored (posixct) idx_calendar ",
         "(cal$anchor is not a Date/POSIXct/yearqtr/yearmon, or cal$posixct ",
         "is FALSE). For non-calendar anchors, use idx_offset_to_pos() with ",
         "a plain numeric offset instead of a date.")
  }
  
  if (!is.null(cal$multi_step)) {
    return(idx_multi_step_date_to_offset(cal, date))
  }
  if (is.na(cal$amount)) {
    # A single compound idx_step repeated every position: strictly
    # monotonic in n, so invertible by binary search rather than walking.
    return(idx_step_date_to_offset(cal, date))
  }
  
  if (is_yearqtr_anchor || is_yearmon_anchor) {
    if (is_yearqtr_anchor && cal$unit != "quarters") {
      stop("idx_to_pos: a yearqtr cal$anchor requires cal$unit = 'quarters'.")
    }
    if (is_yearmon_anchor && cal$unit != "months") {
      stop("idx_to_pos: a yearmon cal$anchor requires cal$unit = 'months'.")
    }
  }
  
  idx_calendar_date_to_pos(cal, date)
}

#' @title Translate a plain numeric offset to an \code{idx_series} integer
#' position
#'
#' @description Inverse of \code{\link{idx_to_date}}'s plain-arithmetic
#' branch (\code{posixct = FALSE}, or a non-Date/POSIXct \code{cal$anchor}):
#' converts a value already expressed in \code{cal$anchor}'s own units/scale
#' (e.g. a picosecond count) into the integer \code{idx_series} position
#' that corresponds to it.
#'
#' Only relevant when you actually have an \code{idx_calendar} with a
#' non-Date/POSIXct \code{anchor} (e.g. a numeric anchor for sub-second or
#' otherwise non-calendar timesteps) and want to look up a position from a
#' value in that anchor's scale. If your series has no calendar meaning at
#' all - no dates, no anchor, nothing to translate - you don't need an
#' \code{idx_calendar} or this function: \code{idx_series} positions
#' (\code{start}/\code{end} arguments throughout this package) are already
#' plain integers, so just use them directly (see \code{\link{idx_range}}).
#'
#' @param cal An \code{idx_calendar} object.
#' @param value A single number, in the same units as \code{cal$anchor}
#' (i.e. \code{cal$anchor + n * cal$amount} for some integer \code{n}), or
#' a vector of such.
#'
#' @returns An integer (vector, matching \code{value}): the \code{idx_series}
#' position(s) corresponding to \code{value}.
#'
#' @examples
#' cal <- idx_calendar(anchor = 0, anchor_pos = 1L, amount = 2.5,
#'                      unit = "picoseconds")
#' idx_offset_to_pos(cal, 12.5)
#'
#' @export
idx_offset_to_pos <- function(cal, value) {
  stopifnot(is_idx_calendar(cal))
  steps <- (value - cal$anchor) / cal$amount
  if (!isTRUE(all.equal(steps, round(steps)))) {
    stop("idx_offset_to_pos: value does not fall on a whole step boundary ",
         "relative to cal$anchor.")
  }
  as.integer(cal$anchor_pos + round(steps))
}

#' @title Compute log growth rate of cumulated dataset
#'
#' @description Helper method to compute the log growth rates of cumulated
#' variables.
#'
#' @param dt Cumulated data series, as an \code{idx_series} with exactly one
#' column.
#' @returns An \code{idx_series} of log growth rates of the cumulated
#' variable inputted via the parameter \code{dt}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1)
#' df2ldl(x)
#'
#' @export
df2ldl <- function(dt) {
  if (!is_idx_series(dt)) {
    stop("dt must be an idx_series object.")
  }
  if (idx_ncol(dt) != 1) {
    stop("dt must only contain 1 data column.")
  }
  lagged <- idx_lag(dt, 1L)
  pos <- idx_positions(dt)
  overlap <- intersect(pos, idx_positions(lagged))
  lag.ov <- lagged[overlap]
  if (any(idx_values(lag.ov) < 0, na.rm = TRUE)) {
    stop("Dataset dt contains negative values.")
  }
  d <- idx_diff(dt, 1L)
  if (any(idx_values(d) < 0, na.rm = TRUE)) {
    stop("Dataset dt has nonpositive increments.")
  }
  lag.aligned <- lagged[idx_positions(d)]
  ldl <- log(idx_values(d) / idx_values(lag.aligned))
  # Leading NA at dt's start position: no preceding value to compute from.
  idx_series(c(NA_real_, ldl), start = dt$start)
}

#' @title Subsetting \code{idx_series} objects given start and end positions
#'
#' @description Helper method to subset an \code{idx_series} for a
#' specified range of integer positions.
#'
#' @param df An \code{idx_series} object, or \code{NULL} if no data are
#' supplied.
#' @param start Start position (integer) of the range.
#' @param end End position (integer) of the range. Defaults to the last
#' position in \code{df}.
#' @returns An \code{idx_series} containing the selected observations, or
#'   \code{NULL} if \code{df} is \code{NULL}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1)
#' get_timeframe(x, 5, 10)
#' get_timeframe(x, 5)
#'
#' @export
get_timeframe <- function(df, start, end = NULL) {
  if (is.null(df)) {
    return(NULL)
  }
  if (!is_idx_series(df)) {
    stop("df is not an idx_series object.")
  }
  rng <- idx_range(df)
  if (is.null(end)) {
    end <- rng[2]
  }
  if (length(start) != 1 || length(end) != 1) {
    stop("start and end must each be a single integer position.")
  }
  start <- max(start, rng[1])
  end <- min(end, rng[2])
  if (start > end) {
    stop("start is after end within the available range of df.")
  }
  df[start:end]
}

#' @title Compute successive increments and log growth rate of 2-variable
#' cumulated dataset
#'
#' @description Helper method to compute the successive increments and log
#' growth rates of cumulated variables. It will compute the successive
#' increments and log cumulative growth rate for each column in the
#' 2-column series, which will then be used to predict or estimate with the
#' leading indicator model.
#'
#' @param data Cumulated data series as an \code{idx_series} with 2 columns:
#' leading indicator and target variable. Can specify which column is
#' leading indicator by \code{LeadIndCol} parameter.
#' @param LeadIndCol Column number of \code{data} that contains the leading
#' indicator. An integer that can only take values 1 (by default) or 2.
#' @returns A list of \code{idx_series} with the original cumulative
#' variables, successive increments and log growth rates: \code{cLead},
#' \code{cTarg}, \code{newLead}, \code{newTarg}, \code{LDLlead},
#' \code{LDLtarg}.
#'
#' @export
add_daily_ldl <- function(data, LeadIndCol = 1) {
  if (!is_idx_series(data)) {
    stop("data is not an idx_series object.")
  }
  if (idx_ncol(data) != 2) {
    stop("Dataset data must contain exactly two series.")
  }
  mat <- idx_values(data)
  if (LeadIndCol == 1) {
    cLead <- idx_series(mat[, 1], start = data$start)
    cTarg <- idx_series(mat[, 2], start = data$start)
  } else if (LeadIndCol == 2) {
    cLead <- idx_series(mat[, 2], start = data$start)
    cTarg <- idx_series(mat[, 1], start = data$start)
  } else {
    stop("LeadIndCol must be an integer, either 1 or 2.")
  }
  
  newLead <- idx_diff(cLead, 1L)
  newTarg <- idx_diff(cTarg, 1L)
  LDLlead <- df2ldl(cLead)
  LDLtarg <- df2ldl(cTarg)
  
  list(
    cLead = cLead, cTarg = cTarg,
    newLead = newLead, newTarg = newTarg,
    LDLlead = LDLlead, LDLtarg = LDLtarg
  )
}

#' @title Reinitialise a series at a given position
#'
#' @description Takes a cumulated series and re-bases it so that it starts
#' from zero at \code{reinit.idx - 1}.
#'
#' @param dt Cumulated data series, as an \code{idx_series} with exactly one
#' column.
#' @param reinit.idx Integer position at which reinitialisation should
#' occur (i.e. \eqn{t=r}, using the notation in the vignette).
#'
#' @returns The reinitialised series, as an \code{idx_series} starting at
#' \code{reinit.idx}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1)
#' reinitialise_dataframe(x, 10)
#'
#' @export
reinitialise_dataframe <- function(dt, reinit.idx) {
  if (!is_idx_series(dt)) {
    stop("dt is not an idx_series object.")
  }
  if (idx_ncol(dt) != 1) {
    stop("dt must only contain 1 data column.")
  }
  rng <- idx_range(dt)
  if (reinit.idx < rng[1] + 1 || reinit.idx > rng[2]) {
    stop("reinit.idx is not present in dt (or has no preceding value).")
  }
  base <- as.numeric(idx_values(dt[reinit.idx - 1]))
  sub <- dt[reinit.idx:rng[2]]
  idx_series(as.numeric(idx_values(sub)) - base, start = reinit.idx)
}

#' @title Return index and value of maximum
#' @description Similar to Python's argmax function.
#' @param x Object to have its maximum found; either an \code{idx_series}
#' or a plain numeric vector.
#' @param decreasing Logical value indicating whether \code{x} should be
#' ordered in decreasing order. Default is \code{TRUE}. Setting this to
#' \code{FALSE} would find the minimum.
#' @returns If \code{x} is an \code{idx_series}, a length-1 \code{idx_series}
#' at the position of the maximum. Otherwise the maximum value.
#' @examples
#' x <- idx_series(cumsum(rpois(30, 5)) + 1, start = 1)
#' argmax(x)
#' @export
argmax <- function(x, decreasing = TRUE) {
  if (is_idx_series(x)) {
    vals <- idx_values(x)
    ord <- order(vals, decreasing = decreasing)[1]
    return(x[idx_positions(x)[ord]])
  }
  return(x[order(x, decreasing = decreasing)[1]])
}

#' @title Compute Mean Absolute Percentage Error (MAPE) for Forecasts Against
#' a Holdout Sample
#'
#' @description This is a helper function that calculates five error metrics
#' of a forecast generated by time series growth curve (tsgc) models. It
#' compares the forecasted values to a holdout sample, providing a measure
#' of forecast accuracy.
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param n.ahead Integer specifying the number of periods to forecast
#' ahead.
#' @param Y An \code{idx_series} object containing the original cumulative
#' dataset.
#'
#' @returns A list containing five error metrics for the forecast, with
#' element names
#' \itemize{
#' \item mape: mean absolute percentage error
#' \item smape: symmetric mean absolute percentage error (between 0 to 100)
#' \item mae: mean absolute error
#' \item rmse: root mean squared error
#' \item coverage: Percentage of holdout sample data points that lie inside
#' the confidence interval for predictions}
#'
#' @export
mapes <- function(res, n.ahead, Y) {
  res$mapes(n.ahead, Y)
}

#' @title Estimate the reproduction number from a fitted model
#'
#' @description Computes the (instantaneous) reproduction number
#' \eqn{R_t = \exp\{g_t \times gen\_ int\}}, where \eqn{g_t} is the
#' filtered or smoothed growth rate of the incidence variable \eqn{y}
#' returned by \code{res$get_gy_ci()}, and \code{gen_int} is the mean
#' generation interval of the disease/process being modelled (in the same
#' integer-position units as \code{res}). This estimate, and its confidence
#' interval, is a deterministic transformation of \code{get_gy_ci()}'s
#' output; \code{res} must already have been produced by \code{estimate()}.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object,
#' obtained from the \code{estimate()} method.
#' @param gen_int The mean generation interval, in the same integer-position
#' units as \code{res} (e.g. if positions are days, this is the mean
#' generation interval in days).
#' @param n.ahead Number of most recent integer positions to return, taken
#' from the end of \code{res}'s filtered/smoothed range. Default is
#' \code{7}.
#' @param smoothed Logical value indicating whether to use the smoothed
#' (\code{TRUE}) or filtered (\code{FALSE}, default) growth rate estimates
#' underlying \eqn{R_t}. Passed through to \code{res$get_gy_ci()}.
#' @param confidence.level Confidence level for the confidence interval
#' around \eqn{R_t}. Default is \code{0.68}, one standard deviation for a
#' normally distributed random variable. Passed through to
#' \code{res$get_gy_ci()}.
#'
#' @returns A data frame containing \code{Position}, \code{fit},
#' \code{lower} and \code{upper}. If \code{res$calendar} is available, a
#' \code{Date} column is included before \code{Position}. Each row gives the
#' estimated reproduction number and its confidence interval over one of the
#' last \code{n.ahead} integer positions.
#'
#' @examples
#' library(tsgc)
#' set.seed(1)
#' Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1)
#' cal <- idx_calendar(anchor = as.Date("2021-01-01"), anchor_pos = 1L,
#'                     amount = 1, unit = "days")
#' model <- SSModelDynamicGompertz$new(Y = Y, q = NULL, end = 100,
#'                                     calendar = cal)
#' res <- estimate(model)
#' estimate_r0(res, gen_int = 5, n.ahead = 7)
#'
#' @export
estimate_r0 <- function(res, gen_int, n.ahead = 7, smoothed = FALSE,
                        confidence.level = 0.68) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  gy.ci <- res$get_gy_ci(smoothed = smoothed, confidence.level = confidence.level)
  rt_mat <- exp(idx_values(gy.ci) * gen_int)
  n <- nrow(rt_mat)
  keep <- seq.int(max(1L, n - n.ahead + 1L), n)
  positions <- gy.ci$start + keep - 1L
  out <- data.frame(
    Position = positions,
    as.data.frame(rt_mat[keep, , drop = FALSE]),
    check.names = FALSE
  )
  if (is_idx_calendar(res$calendar)) {
    out <- data.frame(
      Date = idx_to_date(res$calendar, positions),
      out,
      check.names = FALSE
    )
  }
  out
}

#' @title Walk-Forward Validation for Model Comparison Using Mean Absolute
#' Percentage Error (MAPE)
#'
#' @description This function performs a walk-forward validation to compare
#' forecasting performance across different models specified by the user.
#' It returns a data frame of a user-specified error metric (e.g. MAPE, MAE)
#' for forecasts \code{n.ahead} positions ahead, using the given models with
#' varying end positions.
#'
#' @param Y An \code{idx_series} representing the cumulative data series. If
#' a Leading Indicator model is compared, Y should include columns for both
#' the leading indicator and the target variable. The specific column for
#' the leading indicator can be designated using the \code{LeadIndCol}
#' parameter.
#' @param model_list A list containing \code{SSModelDynamicGompertz} or
#' \code{SSModelLeadingIndicator} objects, to be compared in a cross
#' validation procedure.
#' @param est.end The initial estimation end position for model fitting.
#' Starting from this position, the function re-estimates the model and
#' evaluates the performance for each lag in \code{all_lags} every
#' \code{gap} positions, over \code{n.estimate} steps.
#' @param n.ahead Integer specifying the number of positions to forecast
#' ahead for MAPE evaluation.
#' @param n.estimate Integer indicating the total number of walk-forward
#' validation steps to report.
#' @param gap Integer specifying the position gap between two successive
#' validations, where the model is re-estimated and evaluated during the
#' walk-forward validation.
#' @param xpred_lead.full (Only required for leading indicator models) An
#' \code{idx_series} containing the values of exogenous variables for the
#' leading indicator over the estimation and prediction time frame.
#' @param xpred_targ.full An \code{idx_series} containing the values of
#' exogenous variables for the target variable over the estimation and
#' prediction time frame.
#' @param LeadIndCol (Only required for leading indicator models) Integer
#' representing the column number in \code{Y} that contains the leading
#' indicator.
#' @param criterion A string object indicating how to compare between
#' different models. Available choices are "mape" (by default), "smape",
#' "mae" and "rmse".
#'
#' @returns A table summarizing the chosen error metric for each model in
#' \code{model_list} across the specified positions.
#'
#' @export
cross_val <- function(Y, model_list, est.end, n.ahead = 7, n.estimate = 1, gap = 1,
                      xpred_targ.full = NULL, xpred_lead.full = NULL,
                      LeadIndCol = 1, criterion = "mape") {
  if (!is_idx_series(Y)) {
    stop("Y must be an idx_series object.")
  }
  if (idx_ncol(Y) == 1) {
    Y1 <- Y
  } else if (idx_ncol(Y) == 2) {
    keep_col <- setdiff(1:2, LeadIndCol)
    Y1 <- idx_series(idx_values(Y)[, keep_col], start = Y$start)
  } else {
    stop("Y should not have more than 2 columns.")
  }
  if (length(est.end) != 1 || !isTRUE(all.equal(est.end, as.integer(est.end)))) {
    stop("est.end must be a single integer position.")
  }
  if (n.ahead <= 0) {
    stop("n.ahead must be a positive integer.")
  }
  results <- data.frame(
    Model = names(model_list)
  )
  for (k in 1:n.estimate) {
    index_num <- 1
    for (model in model_list) {
      model$end <- est.end + (k - 1) * gap
      if (inherits(model, "SSModelDynamicGompertz")) {
        model$Y <- get_timeframe(Y1, model$start, model$end)
        if (!is.null(model$xpred)) {
          model$xpred <- get_timeframe(xpred_targ.full, model$start, model$end)
        }
        res <- estimate(model)
        if (res$xpred_logical) {
          res$xpred.new <- xpred_targ.full
        }
        results[index_num, k + 1] <- round(mapes(res, n.ahead, Y1)[[criterion]], 2)
      } else if (inherits(model, "SSModelLeadingIndicator")) {
        if (!is.null(model$xpred_lead)) {
          model$xpred_lead <- xpred_lead.full
        }
        if (!is.null(model$xpred_targ)) {
          model$xpred_targ <- xpred_targ.full
        }
        res <- estimate(model)
        if (res$xpred_logical[1]) {
          res$xpred_lead.new <- xpred_lead.full
        }
        if (res$xpred_logical[2]) {
          res$xpred_targ.new <- xpred_targ.full
        }
        results[index_num, k + 1] <- round(mapes(res, n.ahead, Y)[[criterion]], 2)
      } else {
        stop(paste("Model", index_num, "in model_list is not a SSModelDynamicGompertz or SSModelLeadingIndicator object."))
      }
      index_num <- index_num + 1
    }
  }
  all_ends <- as.character(est.end + c(0:(k - 1)) * gap)
  colnames(results) <- c("Model", all_ends)
  return(results)
}

#' @title Write a selection of relevant results to disk
#'
#' @description Function writes the following results to csv files which get
#' saved in the location specified in \code{res.dir}: forecast new cases or
#' incidence variable, \eqn{y}; the filtered level and slope of \eqn{\ln g},
#' \eqn{\delta} and \eqn{\gamma}; filtered estimates of \eqn{g_y} and the
#' confidence intervals for these estimates.
#'
#' @param res Results object of class \code{FilterResults} or \code{FilterResultsLI},
#' obtained from \samp{estimate()} method.
#' @param res.dir File path to save the results to. A character string.
#' @param n.ahead Number of periods ahead to forecast. A positive integer.
#' @param prefix The prefix to be added to the file names generated. A character string.
#' @param confidence.level Confidence level to use for the confidence interval
#' on the forecasts \eqn{\ln(g_t)}.
#'
#' @importFrom utils write.csv
#' @importFrom stats qnorm
#'
#' @returns A number of csv files saved in the directory specified in
#' \code{res.dir}.
#' @examples
#' # Not run as do not wish to save to local disk when compiling documentation.
#' # Below will run if copied and pasted into console.
#' library(tsgc)
#' library(here)
#'
#' res.dir <- tempdir()
#' data(gauteng,package="tsgc")
#' conv <- xts_to_idx(gauteng)
#' res <- estimate(SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#' calendar = conv$calendar,
#' end = idx_to_pos(conv$calendar, "2020-07-06")))
#'
#' tsgc::write_results(
#' res=res, res.dir = res.dir, prefix="dyn_gompertz",n.ahead = 14,
#' confidence.level = 0.68)
#'
#' @export
write_results <- function(res, res.dir, n.ahead, prefix="", confidence.level=0.68) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")){
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  
  # New cases (delta Y)
  y.hat.diff <- res$predict_level(
    n.ahead = n.ahead,
    confidence.level = confidence.level,
    sea.on = TRUE)
  
  write.csv(
    as.matrix(y.hat.diff),
    row.names = idx_positions(y.hat.diff),
    file = file.path(res.dir, paste(prefix, "cases_fcst.csv", sep="")))
  
  # Filtered slope / level
  y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
  filtered.level <- y.hat.all$level.t.t
  filtered.slope <- y.hat.all$slope.t.t
  a.t.t <- y.hat.all$a.t.t
  P.t.t <- y.hat.all$P.t.t
  idx.slope <- grep("slope", colnames(a.t.t))
  idx.level <- grep("level", colnames(a.t.t))[1]
  gamma.std.err <- sqrt(P.t.t[idx.slope, idx.slope,])
  delta.std.err <- sqrt(P.t.t[idx.level, idx.level,])
  
  gamma <- idx_cbind(filtered.slope, idx_series(gamma.std.err, start = filtered.slope$start))
  delta <- idx_cbind(filtered.level, idx_series(delta.std.err, start = filtered.level$start))
  gamma.mat <- as.matrix(gamma)
  delta.mat <- as.matrix(delta)
  colnames(gamma.mat) <- c("gamma", "std.err")
  colnames(delta.mat) <- c("delta", "std.err")
  
  write.csv(
    gamma.mat,
    row.names = idx_positions(filtered.slope),
    file = file.path(res.dir, paste(prefix, "trend_slope_filt.csv", sep=""))
  )
  write.csv(
    delta.mat,
    row.names = idx_positions(filtered.level),
    file = file.path(res.dir, paste(prefix, "log_gr_level_filt.csv", sep=""))
  )
  
  # Filtered growth rate of new cases (g_y); CI from standard error on the
  # slope component of the state covariance matrix.
  g.y.t.t <- exp(idx_values(filtered.level)) + idx_values(filtered.slope)
  ci <- qnorm((1 - confidence.level) / 2) * gamma.std.err %o% c(1, -1)
  ci_bounds <- as.vector(g.y.t.t) + ci
  gy.ci.mat <- cbind(fit = g.y.t.t, prediction = ci_bounds)
  colnames(gy.ci.mat) <- c("fit", "lower", "upper")
  
  write.csv(
    gy.ci.mat,
    row.names = idx_positions(filtered.level),
    file = file.path(res.dir, paste(prefix, "cases_gr.csv", sep="")))
  
  message("Saved results for: ", substitute(res))
}

#' @title Check whether estimated variance parameters sit at or near a
#' degeneracy boundary
#'
#' @description Flags estimated variance parameters (observation noise
#' \eqn{H}, and any state-innovation variances such as the slope or
#' seasonal components of \eqn{Q}) that are at or below \code{boundary_tol}.
#' A variance parameter this close to zero is a common sign of an
#' unidentified or degenerate parameter in these state-space models.
#'
#' \code{H} can be a scalar (e.g. \code{\link{FilterResults}}, one
#' observation equation) or a \eqn{K \times K} matrix for \eqn{K > 1}
#' (e.g. \code{\link{FilterResultsLI}}'s leading-indicator model, which has
#' two observation equations). This check works identically either way -
#' every element of \code{H} is checked, not just a single assumed-scalar
#' value.
#'
#' @param ... One or more variance parameters (scalars, vectors, or
#' matrices) to check. \code{NA}/\code{NULL} arguments are ignored.
#' @param boundary_tol Boundary threshold. Defaults to \code{1e-6}.
#'
#' @returns Logical scalar: \code{TRUE} if any supplied variance parameter
#' is non-missing and at or below \code{boundary_tol}, \code{FALSE}
#' otherwise (including when every supplied parameter is entirely
#' missing/\code{NULL}).
#' @export
check_variance_boundary <- function(..., boundary_tol = 1e-6) {
  params <- list(...)
  any(vapply(params, function(p) {
    if (is.null(p)) return(FALSE)
    p <- as.numeric(p)
    if (anyNA(p)) p <- p[!is.na(p)]
    if (length(p) == 0) return(FALSE)
    any(p <= boundary_tol)
  }, logical(1)))
}


#' @title Calling print method for classes in tsgc
#'
#' @description Accessor method to print a short description for the objects of
#' `SSModelLeadingIndicator` class
#'
#' @param x A `SSModelLeadingIndicator` object
#' @param ... Additional arguments.
#' 
#' @method print SSModelLeadingIndicator
#' 
#' @examples
#' library(tsgc)
#' data(england, package = "tsgc")
#' conv <- xts_to_idx(england[, 1:2])
#'
#' # Specify a model
#' out_eng <- SSModelLeadingIndicator(
#'   Y = conv$series, n.lag = 4, sea.period = 7, LeadIndCol = 1,
#'   calendar = conv$calendar,
#'   start = idx_to_pos(conv$calendar, "2021-04-30"),
#'   end = idx_to_pos(conv$calendar, "2021-07-24"))
#' 
#' # Print a short description of the model object
#' print(out_eng)
#' 
#' 
#' @export
print.SSModelLeadingIndicator <- function(x, ...) {
  x$print()
}


#' @title Calling summary method for classes in tsgc
#'
#' @description Accessor method to show a summary for the objects of
#' `SSModelLeadingIndicator` class
#'
#' @param object A `SSModelLeadingIndicator` object
#' @param ... Additional arguments.
#' @method summary SSModelLeadingIndicator
#' 
#' @examples
#' library(tsgc)
#' data(england, package = "tsgc")
#' conv <- xts_to_idx(england[, 1:2])
#'
#' # Specify a model
#' out_eng <- SSModelLeadingIndicator(
#'   Y = conv$series, n.lag = 4, sea.period = 7, LeadIndCol = 1,
#'   calendar = conv$calendar,
#'   start = idx_to_pos(conv$calendar, "2021-04-30"),
#'   end = idx_to_pos(conv$calendar, "2021-07-24"))
#' 
#' summary(out_eng)
#' 
#' @export
summary.SSModelLeadingIndicator <- function(object, ...) {
  object$summary()
}

#' @title Calling print method for SSModelDynamicGompertz class
#'
#' @description Accessor method to print a short description for the objects of
#' `SSModelDynamicGompertz` class
#'
#' @param x A `SSModelDynamicGompertz` object
#' @param ... Additional arguments.
#' @method print SSModelDynamicGompertz
#' 
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-06"))
#' 
#' # Print a short description of the model object
#' print(model)
#' 
#' @export
print.SSModelDynamicGompertz <- function(x, ...) {
  x$print()
}

#' @title Calling summary method for SSModelDynamicGompertz class
#'
#' @description Accessor method to show a summary for the objects of
#' `SSModelDynamicGompertz` class
#'
#' @param object A `SSModelDynamicGompertz` object
#' @param ... Additional arguments.
#' @method summary SSModelDynamicGompertz
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-06"))
#' 
#' # Show summary of the model object
#' summary(model)
#' 
#' @export
summary.SSModelDynamicGompertz <- function(object, ...) {
  object$summary()
}

#' @title Calling summary method for FilterResults
#'
#' @description Accessor method to show a summary for the objects of
#' `FilterResults` class
#'
#' @param object A `FilterResults` object
#' @param ... Additional arguments.
#' @method summary FilterResults
#' 
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Return KFS object in output of res
#' summary(res)
#' 
#' @export
summary.FilterResults <- function(object, ...) {
  object$summary()
}

#' @title Calling print method for FilterResults class
#'
#' @description Accessor method to print a short description for the objects of
#' `FilterResults` class
#'
#' @param x A `FilterResults` object
#' @param ... Additional arguments.
#' @method print FilterResults
#' 
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Return short description of fitted model
#' print(res)
#' 
#' @export
print.FilterResults <- function(x, ...) {
  x$print()
}

#' @title Calling summary method for FilterResultsLI
#'
#' @description Accessor method to show a summary for the objects of
#' `FilterResultsLI` class
#'
#' @param object A `FilterResultsLI` object
#' @param ... Additional arguments.
#' @method summary FilterResultsLI
#' 
#' @examples
#' library(tsgc)
#' data(england, package = "tsgc")
#' conv <- xts_to_idx(england[, 1:2])
#'
#' out_eng <- SSModelLeadingIndicator(
#'   Y = conv$series, n.lag = 4, sea.period = 7, LeadIndCol = 1,
#'   calendar = conv$calendar,
#'   start = idx_to_pos(conv$calendar, "2021-04-30"),
#'   end = idx_to_pos(conv$calendar, "2021-07-24"))
#' 
#' res_eng<-estimate(out_eng)
#' summary(res_eng)
#' 
#' @export
summary.FilterResultsLI <- function(object, ...) {
  object$summary()
}

#' @title Calling print method for FilterResultsLI class
#'
#' @description Accessor method to print a short description for the objects of
#' `FilterResultsLI` class
#'
#' @param x A `FilterResultsLI` object
#' @param ... Additional arguments.
#' @method print FilterResultsLI
#' 
#' @examples
#' library(tsgc)
#' data(england, package = "tsgc")
#' conv <- xts_to_idx(england[, 1:2])
#'
#' out_eng <- SSModelLeadingIndicator(
#'   Y = conv$series, n.lag = 4, sea.period = 7, LeadIndCol = 1,
#'   calendar = conv$calendar,
#'   start = idx_to_pos(conv$calendar, "2021-04-30"),
#'   end = idx_to_pos(conv$calendar, "2021-07-24"))
#' 
#' res_eng<-estimate(out_eng)
#' print(res_eng) 
#' 
#' @export
print.FilterResultsLI <- function(x, ...) {
  x$print()
}



#' @title Report fitted diagnostics for a FilterResults/FilterResultsLI object
#'
#' @description \code{summary()} on \code{\link{FilterResults}}/
#' \code{\link{FilterResultsLI}} reports estimated variance parameters and
#' model states, but does not report the model's log-likelihood, a
#' boundary/degeneracy check across all estimated variance parameters, or
#' residual diagnostics. \code{print_model_diagnostics()} reports these
#' three additional, complementary diagnostics; it does not replace or
#' duplicate what \code{summary()} already prints.
#'
#' @param res A \code{FilterResults} or \code{FilterResultsLI} object.
#' @param boundary_tol Boundary threshold passed to
#' \code{\link{check_variance_boundary}}. Defaults to \code{1e-6}.
#'
#' @returns \code{NULL}, invisibly. Called for its printed output.
#' @importFrom stats residuals sd
#' @export
print_model_diagnostics <- function(res, boundary_tol = 1e-6) {
  if (!inherits(res, "FilterResults") && !inherits(res, "FilterResultsLI")) {
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  cat("---- Model diagnostics ----\n")
  
  # Parameter estimates.
  tryCatch(res$print_estimation_results(), error = function(e) {
    cat("  (print_estimation_results unavailable: ", conditionMessage(e), ")\n", sep = "")
  })
  
  # Log-likelihood of the fitted state-space model, when available.
  ll <- tryCatch(stats::logLik(res$output$model), error = function(e) NA)
  cat("  Log-likelihood:", if (is.na(ll)) "unavailable" else format(ll, digits = 6), "\n")
  
  # Flag estimated variances at or near a degeneracy boundary. H may be a
  # K x K matrix for multi-equation models (e.g. FilterResultsLI).
  H  <- tryCatch(res$output$model$H[, , 1], error = function(e) NULL)
  Qg <- tryCatch(res$output$model$Q[2, 2, 1], error = function(e) NULL)
  boundary_hit <- check_variance_boundary(H, Qg, boundary_tol = boundary_tol)
  cat("  Boundary check (variance <= ", boundary_tol, "): ",
      if (boundary_hit) "FLAGGED - one or more variances at/near zero" else "ok",
      "\n", sep = "")
  
  # Residual/innovation diagnostics from the KFS object, when present.
  resid <- tryCatch(residuals(res$output, type = "recursive"), error = function(e) NULL)
  if (!is.null(resid)) {
    resid <- resid[is.finite(resid)]
    cat("  Recursive residuals: mean =", format(mean(resid), digits = 4),
        ", sd =", format(sd(resid), digits = 4), "\n")
  } else {
    cat("  Recursive residuals: unavailable\n")
  }
  cat("----------------------------\n")
  invisible(NULL)
}