#' @title Canonical form of a calendar unit string
#'
#' @description Internal helper. Recognises the eight standard calendar
#' units - seconds, minutes, hours, days, weeks, months, quarters, years
#' (singular or plural, case-insensitive) - and returns their canonical
#' singular and plural forms. Returns \code{NULL} for any other unit
#' string. This is the single source of truth for unit-name recognition,
#' used both to build an \code{\link{idx_step}} from a plain
#' \code{amount}/\code{unit} pair (see \code{\link{idx_calendar}}) and to
#' map a unit to the \code{by} string expected by
#' \code{seq.Date}/\code{seq.POSIXt} (see \code{\link{idx_calendar_by_unit}}).
#'
#' @param unit A single string.
#' @returns \code{NULL} if \code{unit} is not one of the eight standard
#' units, otherwise a list with elements \code{plural} (matching an
#' \code{idx_step} field name, e.g. \code{"seconds"}) and \code{singular}
#' (matching a \code{seq.Date}/\code{seq.POSIXt} \code{by} unit, e.g.
#' \code{"sec"}).
#' @keywords internal
#' @noRd
idx_canonical_unit <- function(unit) {
  u <- tolower(unit)
  switch(u,
         second = , seconds = list(plural = "seconds", singular = "sec"),
         minute = , minutes = list(plural = "minutes", singular = "min"),
         hour = , hours = list(plural = "hours", singular = "hour"),
         day = , days = list(plural = "days", singular = "day"),
         week = , weeks = list(plural = "weeks", singular = "week"),
         month = , months = list(plural = "months", singular = "month"),
         quarter = , quarters = list(plural = "quarters", singular = "quarter"),
         year = , years = list(plural = "years", singular = "year"),
         NULL
  )
}

#' @title Classify an anchor's calendar kind
#'
#' @description Internal helper. Determines which of the four supported
#' anchor kinds (\code{yearqtr}, \code{yearmon}, \code{Date}/\code{POSIXct},
#' or none of the above - e.g. a plain number) \code{anchor} is, in one
#' place. Used everywhere an \code{idx_step}/\code{idx_calendar} function
#' needs to branch on anchor kind, instead of each site re-deriving the
#' same set of \code{inherits()} checks independently.
#'
#' @param anchor An anchor value (see \code{\link{idx_calendar}}).
#' @returns A list with logical elements \code{is_yearqtr},
#' \code{is_yearmon}, \code{is_posixct} (\code{TRUE} for \code{POSIXct}
#' only), \code{is_date_posix} (\code{TRUE} for \code{Date} or
#' \code{POSIXct}), and \code{is_calendar} (\code{TRUE} if any of
#' \code{is_yearqtr}/\code{is_yearmon}/\code{is_date_posix} is
#' \code{TRUE}).
#' @keywords internal
#' @noRd
idx_anchor_kind <- function(anchor) {
  is_yearqtr <- inherits(anchor, "yearqtr")
  is_yearmon <- inherits(anchor, "yearmon")
  is_posixct <- inherits(anchor, "POSIXct")
  is_date_posix <- inherits(anchor, "Date") || is_posixct
  list(is_yearqtr = is_yearqtr, is_yearmon = is_yearmon,
       is_posixct = is_posixct, is_date_posix = is_date_posix,
       is_calendar = is_yearqtr || is_yearmon || is_date_posix)
}

#' @title A compound calendar step
#'
#' @description \code{idx_step} describes the size of a single step between
#' consecutive \code{idx_series} integer positions as a combination of
#' calendar-relative and fixed-duration components: years, quarters,
#' months, weeks, days, hours, minutes, and seconds. Each component is
#' independent - they are not normalized or collapsed into one another
#' (e.g. \code{idx_step(months = 1, days = 31)} is not simplified to two
#' months; it is applied as "one month, then 31 days"). This lets a step
#' express things a single \code{amount}/\code{unit} pair cannot, such as
#' "3 seconds after each quarter" (\code{idx_step(quarters = 1, seconds = 3)})
#' or "a month and 13 seconds" (\code{idx_step(months = 1, seconds = 13)}).
#'
#' When applied (see \code{\link{idx_step_add}}), components are added to
#' an anchor date in a fixed order, from largest calendar-relative unit to
#' smallest fixed-duration unit: years, then quarters, then months (each
#' calendar-aware, i.e. respecting variable month/quarter length), then
#' weeks, days, hours, minutes, seconds (each fixed-duration). This order
#' is applied once per full step multiple - i.e. for \code{n} steps from
#' the anchor, every component is scaled by \code{n} before being added,
#' not added once and then repeated.
#'
#' A step with any of \code{hours}, \code{minutes}, or \code{seconds}
#' non-zero requires a \code{POSIXct} anchor (not a plain \code{Date} or
#' \code{yearqtr}/\code{yearmon}), since those units need a time-of-day
#' component to be meaningful.
#'
#' @param years,quarters,months,weeks,days,hours,minutes,seconds Single
#' non-negative, finite numbers giving the count of each component in one
#' step. All default to \code{0}. At least one must be non-zero.
#'
#' @returns An object of class \code{idx_step}.
#'
#' @examples
#' # A single quarter, nudged 3 seconds later - e.g. data recorded 3
#' # seconds after each quarter boundary.
#' idx_step(quarters = 1, seconds = 3)
#'
#' # A month and 13 seconds - the extra seconds accumulate with each step,
#' # same as the calendar-relative part.
#' idx_step(months = 1, seconds = 13)
#'
#' # Equivalent to amount = 1, unit = "days".
#' idx_step(days = 1)
#'
#' @export
idx_step <- function(years = 0, quarters = 0, months = 0, weeks = 0,
                     days = 0, hours = 0, minutes = 0, seconds = 0) {
  fields <- list(years = years, quarters = quarters, months = months,
                 weeks = weeks, days = days, hours = hours,
                 minutes = minutes, seconds = seconds)
  for (nm in names(fields)) {
    v <- fields[[nm]]
    if (length(v) != 1 || !is.numeric(v) || !is.finite(v) || v < 0) {
      stop(nm, " must be a single non-negative, finite number.")
    }
  }
  if (all(vapply(fields, function(v) isTRUE(all.equal(v, 0)), logical(1)))) {
    stop("idx_step: at least one of years, quarters, months, weeks, days, ",
         "hours, minutes, seconds must be non-zero.")
  }
  structure(fields, class = "idx_step")
}

methods::setOldClass("idx_step")

#' @title Test whether an object is an \code{idx_step}
#' @param x Object to test.
#' @returns Logical.
#' @export
is_idx_step <- function(x) {
  inherits(x, "idx_step")
}

#' @title Test whether an \code{idx_step} has any sub-day component
#' @description Internal helper. Returns \code{TRUE} if any of
#' \code{hours}/\code{minutes}/\code{seconds} is non-zero, meaning the step
#' requires a \code{POSIXct} (not plain \code{Date}) anchor to apply.
#' @param step An \code{idx_step} object.
#' @returns Logical.
#' @keywords internal
#' @noRd
idx_step_needs_posixct <- function(step) {
  stopifnot(is_idx_step(step))
  any(vapply(c("hours", "minutes", "seconds"),
             function(nm) !isTRUE(all.equal(step[[nm]], 0)), logical(1)))
}

#' @title Human-readable label for an \code{idx_step}
#'
#' @description Internal helper. Builds the "N unit, N unit, ..." label
#' used both by \code{\link{print.idx_step}} and by the plotting layer's
#' calendar info box (see \code{idx_info_box_caption()} in plotting.R),
#' listing only the non-zero components of \code{step}.
#'
#' @param step An \code{idx_step} object.
#' @returns A single string, e.g. \code{"1 quarters, 3 seconds"}.
#' @keywords internal
#' @noRd
idx_step_label <- function(step) {
  stopifnot(is_idx_step(step))
  parts <- vapply(names(step), function(nm) {
    if (!isTRUE(all.equal(step[[nm]], 0))) paste(step[[nm]], nm) else NA_character_
  }, character(1))
  paste(parts[!is.na(parts)], collapse = ", ")
}

#' @title Print an \code{idx_step}
#' @param x An \code{idx_step} object.
#' @param ... Unused.
#' @export
print.idx_step <- function(x, ...) {
  cat("<idx_step> ", idx_step_label(x), "\n", sep = "")
  invisible(x)
}

#' @title Apply an \code{idx_step} to an anchor date, scaled by a step count
#'
#' @description Adds \code{n} whole steps of \code{step} to \code{anchor},
#' one component at a time, from largest calendar-relative unit to
#' smallest fixed-duration unit (years, quarters, months, weeks, days,
#' hours, minutes, seconds), each scaled by \code{n} before being added.
#' Calendar-relative components (years/quarters/months) use calendar-aware
#' stepping (\code{\link{seq.Date}}/\code{\link{seq.POSIXt}}, or the
#' equivalent fractional-year arithmetic for \code{yearqtr}/\code{yearmon}
#' anchors); fixed-duration components (weeks/days/hours/minutes/seconds)
#' use plain duration addition.
#'
#' @param anchor A \code{Date}, \code{POSIXct}, \code{yearqtr}, or
#' \code{yearmon} to step from.
#' @param step An \code{idx_step} object.
#' @param n An integer (or vector of integers): the number of whole steps
#' to apply. May be negative.
#' @returns A vector the same length as \code{n}, the same class as
#' \code{anchor} (or \code{POSIXct}, if a sub-day component forces
#' promotion from \code{Date}).
#' @export
idx_step_add <- function(anchor, step, n) {
  stopifnot(is_idx_step(step))
  n <- as.integer(round(n))
  
  kind <- idx_anchor_kind(anchor)
  is_yearqtr <- kind$is_yearqtr
  is_yearmon <- kind$is_yearmon
  is_date_posix <- kind$is_date_posix
  
  if (idx_step_needs_posixct(step)) {
    if (is_yearqtr || is_yearmon) {
      stop("idx_step_add: a step with hours/minutes/seconds requires a ",
           "Date or POSIXct anchor, not yearqtr/yearmon (no time-of-day).")
    }
    if (!is_date_posix) {
      stop("idx_step_add: a step with hours/minutes/seconds requires a ",
           "Date or POSIXct anchor.")
    }
    anchor <- as.POSIXct(anchor)
  }
  if ((!isTRUE(all.equal(step$quarters, 0)) || !isTRUE(all.equal(step$months, 0))) &&
      !is_yearqtr && !is_yearmon && !is_date_posix) {
    stop("idx_step_add: quarters/months components require a Date, ",
         "POSIXct, yearqtr, or yearmon anchor.")
  }
  if (is_yearqtr && !isTRUE(all.equal(step$months, 0))) {
    stop("idx_step_add: a 'months' component is not supported for a ",
         "yearqtr anchor (use 'quarters' instead).")
  }
  if (is_yearmon && !isTRUE(all.equal(step$quarters, 0))) {
    stop("idx_step_add: a 'quarters' component is not supported for a ",
         "yearmon anchor (use 'months' instead).")
  }
  if ((is_yearqtr || is_yearmon) &&
      any(vapply(c("weeks", "days", "hours", "minutes", "seconds"),
                 function(nm) !isTRUE(all.equal(step[[nm]], 0)), logical(1)))) {
    stop("idx_step_add: weeks/days/hours/minutes/seconds components are ",
         "not supported for a yearqtr/yearmon anchor (no sub-month resolution).")
  }
  
  step_one <- function(a, ni) {
    if (is_yearqtr || is_yearmon) {
      if (!isTRUE(all.equal(step$years, 0))) a <- a + ni * step$years
      if (!isTRUE(all.equal(step$quarters, 0))) a <- a + ni * step$quarters * 0.25
      if (!isTRUE(all.equal(step$months, 0))) a <- a + ni * step$months * (1/12)
      return(a)
    }
    if (!isTRUE(all.equal(step$years, 0))) {
      a <- seq(a, by = paste(ni * step$years, "year"), length.out = 2)[2]
    }
    if (!isTRUE(all.equal(step$quarters, 0))) {
      a <- seq(a, by = paste(ni * step$quarters, "quarter"), length.out = 2)[2]
    }
    if (!isTRUE(all.equal(step$months, 0))) {
      a <- seq(a, by = paste(ni * step$months, "month"), length.out = 2)[2]
    }
    fixed_days <- ni * step$weeks * 7 + ni * step$days
    fixed_secs <- ni * step$hours * 3600 + ni * step$minutes * 60 + ni * step$seconds
    if (inherits(a, "POSIXct")) {
      a <- a + fixed_days * 86400 + fixed_secs
    } else {
      # fixed_secs is always 0 here for a plain Date anchor.
      a <- a + fixed_days
    }
    a
  }
  
  out <- vapply(n, function(ni) as.numeric(step_one(anchor, ni)), FUN.VALUE = numeric(1))
  class(out) <- class(anchor)
  attr(out, "tzone") <- attr(anchor, "tzone")
  out
}


#' @title A repeating cycle of distinct \code{idx_step}s
#'
#' @description \code{multi_step_pattern} describes a repeating cycle
#' where each slot can be a completely different \code{\link{idx_step}} -
#' not just a different multiple of one fixed step. This is for cases the
#' numeric \code{pattern} argument of \code{\link{idx_calendar}}/
#' \code{\link{idx_calendar_step}} cannot express, because that mechanism
#' only ever scales a single \code{amount}/\code{unit} (or a single
#' compound \code{idx_step}) by different multipliers around the cycle -
#' it cannot mix step *kinds*. For example, "3 days, 3 days, then 1 month,
#' repeating" needs two of the three slots to be fixed-duration and one to
#' be calendar-relative; \code{multi_step_pattern} lets each slot carry
#' its own \code{idx_step}:
#' \code{multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))}.
#'
#' Use this via \code{\link{idx_calendar_multi_step}} (a third,
#' most-general constructor for \code{\link{idx_calendar}}, alongside the
#' primary \code{amount}/\code{unit} constructor and the
#' \code{\link{idx_calendar_step}} compound-step constructor). The
#' existing numeric \code{pattern} mechanism is unchanged and remains the
#' simpler choice whenever every slot in the cycle is just a different
#' multiple of the same single step.
#'
#' @param ... One or more \code{\link{idx_step}} objects, given in the
#' order they repeat.
#'
#' @returns An object of class \code{multi_step_pattern}: a list of
#' \code{idx_step} objects.
#'
#' @examples
#' # 3 days, 3 days, then 1 month - repeating.
#' multi_step_pattern(idx_step(days = 3), idx_step(days = 3), idx_step(months = 1))
#'
#' @export
multi_step_pattern <- function(...) {
  steps <- list(...)
  if (length(steps) < 1) {
    stop("multi_step_pattern: at least one idx_step is required.")
  }
  if (!all(vapply(steps, is_idx_step, logical(1)))) {
    stop("multi_step_pattern: all arguments must be idx_step objects; see idx_step().")
  }
  structure(steps, class = "multi_step_pattern")
}

methods::setOldClass("multi_step_pattern")

#' @title Test whether an object is a \code{multi_step_pattern}
#' @param x Object to test.
#' @returns Logical.
#' @export
is_multi_step_pattern <- function(x) {
  inherits(x, "multi_step_pattern")
}

#' @title Print a \code{multi_step_pattern}
#' @param x A \code{multi_step_pattern} object.
#' @param ... Unused.
#' @export
print.multi_step_pattern <- function(x, ...) {
  cat("<multi_step_pattern> (", length(x), " slots)\n", sep = "")
  for (i in seq_along(x)) {
    cat("  [", i, "] ", sep = "")
    print(x[[i]])
  }
  invisible(x)
}


#' @title Validate fields common to all \code{idx_calendar} constructors
#'
#' @description Internal helper. Validates \code{anchor_name},
#' \code{anchor_pos}, \code{pattern_start} (against a supplied cycle
#' length), and \code{posixct} - the four fields checked identically by
#' \code{\link{idx_calendar}}, \code{\link{idx_calendar_step}}, and
#' \code{\link{idx_calendar_multi_step}}, regardless of which of
#' \code{amount}/\code{unit}, \code{step}, or \code{multi_step} the
#' constructor otherwise takes. Stops with an informative error if any
#' check fails; returns \code{NULL} (invisibly) otherwise.
#'
#' @param anchor_name,anchor_pos,pattern_start,posixct As in
#' \code{\link{idx_calendar}}.
#' @param cycle_len A single positive integer: the length of the cycle
#' \code{pattern_start} must index into (\code{length(pattern)} for
#' \code{\link{idx_calendar}}/\code{\link{idx_calendar_step}}, or
#' \code{length(multi_step)} for \code{\link{idx_calendar_multi_step}}).
#' @keywords internal
#' @noRd
idx_calendar_validate_common <- function(anchor_name, anchor_pos, pattern_start,
                                         posixct, cycle_len) {
  if (!is.null(anchor_name) && (length(anchor_name) != 1 || !is.character(anchor_name))) {
    stop("anchor_name must be NULL or a single string.")
  }
  if (length(anchor_pos) != 1 || !isTRUE(all.equal(anchor_pos, as.integer(anchor_pos)))) {
    stop("anchor_pos must be a single integer.")
  }
  if (length(pattern_start) != 1 || !isTRUE(all.equal(pattern_start, as.integer(pattern_start))) ||
      pattern_start < 1 || pattern_start > cycle_len) {
    stop("pattern_start must be a single integer between 1 and ", cycle_len, ".")
  }
  if (length(posixct) != 1 || !is.logical(posixct) || is.na(posixct)) {
    stop("posixct must be a single logical (TRUE/FALSE).")
  }
  invisible(NULL)
}

#' @title Construct an \code{idx_calendar}
#'
#' @description \code{idx_calendar} describes how to translate the integer
#' positions used by \code{\link{idx_series}} back into calendar time. It
#' holds no data of its own and is not used by the statistical engine of
#' this package (estimation, filtering, forecasting); it exists for use at
#' the edges of the package - plotting, printing, reporting - where a human
#' needs to see a date rather than an integer position.
#'
#' An \code{idx_calendar} is defined by:
#' \itemize{
#'   \item \code{anchor}: a reference point in calendar time (e.g. a
#'   \code{Date}, or a plain number if the series has no true calendar
#'   meaning - e.g. picoseconds since some experiment start).
#'   \item \code{anchor_pos}: the integer position (in the same units as
#'   \code{idx_series$start}) that lines up with \code{anchor}.
#'   \item \code{amount}/\code{unit}: the size of a single step between
#'   consecutive integer positions, e.g. \code{amount = 2.5, unit = "hours"}
#'   or \code{amount = 1, unit = "days"}. \code{amount} may be non-integer.
#'   \item \code{pattern}: a repeating vector of step multipliers, used
#'   when consecutive integer positions do not correspond to evenly spaced
#'   calendar steps. For example, business days (weekends skipped) would
#'   be \code{pattern = c(1, 1, 1, 1, 3)}: four single-day steps, Mon-Thu,
#'   then a 3-day step from Fri to the following Mon. Defaults to
#'   \code{1}, i.e. every step is a single, uniform \code{amount}.
#'   \item \code{pattern_start}: the index into \code{pattern} that
#'   \code{anchor_pos} sits at, i.e. which position in the repeating cycle
#'   the anchor falls on. Defaults to \code{1}.
#' }
#'
#' This object deliberately does no calendar arithmetic itself beyond plain
#' numeric offsets (number of steps x \code{amount}, in \code{unit}s). Any
#' interpretation of that numeric offset as a true calendar date/time
#' (including POSIXct-specific concerns such as DST or variable month
#' lengths) is left to the plotting/translation layer, which is expected to
#' use \code{idx_calendar_offset()} together with the \code{posixct} flag on
#' this object to decide how to turn the offset into an actual date.
#'
#' @param anchor A single reference point in calendar time. Typically a
#' \code{Date} or \code{POSIXct}, but may be a plain number for series with
#' no true calendar meaning.
#' @param anchor_pos A single integer: the \code{idx_series} position that
#' \code{anchor} corresponds to.
#' @param amount A single positive number: the size of one step between
#' consecutive integer positions, in \code{unit}s. May be non-integer (e.g.
#' \code{2.5} for "two and a half hours").
#' @param unit A single string naming the unit that \code{amount} is
#' measured in (e.g. \code{"seconds"}, \code{"hours"}, \code{"days"},
#' \code{"weeks"}, \code{"months"}, \code{"quarters"}, \code{"years"}, or
#' any other unit meaningful to the caller/plotting layer). This package
#' places no constraints on the set of valid units; interpretation is left
#' entirely to the consumer (e.g. the plotting layer's POSIXct conversion).
#' @param pattern A numeric vector of step multipliers, recycled
#' cyclically across consecutive integer positions. Defaults to \code{1}
#' (every step is a uniform \code{amount}). For example, business days
#' (weekends skipped) would use \code{c(1, 1, 1, 1, 3)}.
#' @param pattern_start A single positive integer: the index into
#' \code{pattern} at which \code{anchor_pos} falls. Defaults to \code{1}.
#' @param posixct A single logical. If \code{TRUE}, hints to the plotting/
#' translation layer that \code{anchor} should be treated as (or coerced
#' to) a \code{POSIXct}/\code{Date}, and that calendar-aware arithmetic
#' (e.g. \code{seq.POSIXt}, month/quarter-aware stepping) should be used to
#' interpret offsets rather than plain numeric addition. Defaults to
#' \code{FALSE}. This object does not act on the flag itself; it only
#' stores it for the translation layer to consult.
#'
#' @returns An object of class \code{idx_calendar}.
#'
#' @examples
#' # Daily series, anchored so that position 1 is 2024-01-01.
#' cal <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
#'                      amount = 1, unit = "days", posixct = TRUE)
#' idx_to_date(cal, 1:5)
#'
#' # Business-day series (weekends skipped), anchor at a Monday.
#' cal_bd <- idx_calendar(anchor = as.Date("2024-01-01"), anchor_pos = 1L,
#'                         amount = 1, unit = "days",
#'                         pattern = c(1, 1, 1, 1, 3), pattern_start = 1L,
#'                         posixct = TRUE)
#' idx_to_date(cal_bd, 1:10)
#'
#' # Non-calendar series: picosecond timesteps, arbitrary numeric anchor.
#' cal_ps <- idx_calendar(anchor = 0, anchor_pos = 1L,
#'                         amount = 2.5, unit = "picoseconds")
#' idx_to_date(cal_ps, 1:5)
#'
#' @param anchor_name An optional single string giving a human-readable
#' name for \code{anchor} (e.g. \code{"outbreak start"}, \code{"launch"}).
#' When supplied, the plotting layer can use it to label axes/info boxes
#' as e.g. "steps from launch" or "time since launch" instead of the
#' generic "steps from anchor". Defaults to \code{NULL} (no name; generic
#' wording is used).
#'
#' @export
idx_calendar <- function(anchor, anchor_pos = 1L, amount = 1, unit = "days",
                         pattern = 1, pattern_start = 1L, posixct = FALSE,
                         anchor_name = NULL) {
  if (length(amount) != 1 || !is.numeric(amount) || !is.finite(amount) || amount <= 0) {
    stop("amount must be a single positive, finite number.")
  }
  if (length(unit) != 1 || !is.character(unit)) {
    stop("unit must be a single string.")
  }
  if (!is.numeric(pattern) || length(pattern) < 1 || any(!is.finite(pattern)) || any(pattern <= 0)) {
    stop("pattern must be a numeric vector of one or more positive, finite values.")
  }
  idx_calendar_validate_common(anchor_name, anchor_pos, pattern_start, posixct,
                               cycle_len = length(pattern))
  canonical <- idx_canonical_unit(unit)
  step_obj <- if (!is.null(canonical)) {
    args <- setNames(list(amount), canonical$plural)
    do.call(idx_step, args)
  } else {
    NULL
  }
  structure(
    list(
      anchor = anchor,
      anchor_pos = as.integer(anchor_pos),
      amount = amount,
      unit = unit,
      step = step_obj,
      multi_step = NULL,
      pattern = as.numeric(pattern),
      pattern_start = as.integer(pattern_start),
      posixct = posixct,
      anchor_name = anchor_name
    ),
    class = "idx_calendar"
  )
}

#' @title Construct an \code{idx_calendar} from a compound \code{idx_step}
#'
#' @description Secondary constructor for \code{\link{idx_calendar}},
#' taking a compound \code{\link{idx_step}} (a combination of years,
#' quarters, months, weeks, days, hours, minutes, seconds) as the step
#' size instead of a single \code{amount}/\code{unit} pair. Use this when
#' a single unit cannot express the step - e.g. "3 seconds after each
#' quarter" or "a month and 13 seconds". For a plain single-unit step, the
#' primary \code{\link{idx_calendar}} constructor is simpler.
#'
#' @param anchor A single reference point in calendar time (see
#' \code{\link{idx_calendar}}).
#' @param step An \code{\link{idx_step}} object giving the size of one
#' step between consecutive integer positions.
#' @param anchor_pos,pattern,pattern_start,posixct,anchor_name As in
#' \code{\link{idx_calendar}}.
#'
#' @returns An object of class \code{idx_calendar}.
#'
#' @examples
#' # Data recorded 3 seconds after each quarter boundary.
#' cal <- idx_calendar_step(
#'   anchor = as.POSIXct("2024-01-01 00:00:03", tz = "UTC"),
#'   step = idx_step(quarters = 1, seconds = 3),
#'   posixct = TRUE
#' )
#' idx_to_date(cal, 1:4)
#'
#' @export
idx_calendar_step <- function(anchor, step, anchor_pos = 1L, pattern = 1,
                              pattern_start = 1L, posixct = FALSE,
                              anchor_name = NULL) {
  if (!is_idx_step(step)) {
    stop("step must be an idx_step object; see idx_step().")
  }
  if (!is.numeric(pattern) || length(pattern) < 1 || any(!is.finite(pattern)) || any(pattern <= 0)) {
    stop("pattern must be a numeric vector of one or more positive, finite values.")
  }
  idx_calendar_validate_common(anchor_name, anchor_pos, pattern_start, posixct,
                               cycle_len = length(pattern))
  structure(
    list(
      anchor = anchor,
      anchor_pos = as.integer(anchor_pos),
      amount = NA_real_,
      unit = NA_character_,
      step = step,
      multi_step = NULL,
      pattern = as.numeric(pattern),
      pattern_start = as.integer(pattern_start),
      posixct = posixct,
      anchor_name = anchor_name
    ),
    class = "idx_calendar"
  )
}

methods::setOldClass("idx_calendar")

#' @title Construct an \code{idx_calendar} from a \code{multi_step_pattern}
#'
#' @description Third, most general constructor for
#' \code{\link{idx_calendar}}. Use this when consecutive positions cycle
#' through genuinely different kinds of steps - not just different
#' multiples of one step (that case is covered by the numeric
#' \code{pattern} argument of \code{\link{idx_calendar}}/
#' \code{\link{idx_calendar_step}}). For example, "3 days, 3 days, then 1
#' month, repeating" needs a \code{\link{multi_step_pattern}} since two
#' slots are fixed-duration and one is calendar-relative.
#'
#' Unlike the numeric \code{pattern} mechanism (which pre-computes a
#' closed-form pattern-weighted step count, see
#' \code{\link{idx_calendar_offset}}), a \code{multi_step_pattern} is
#' applied by walking one slot at a time from \code{anchor_pos} to the
#' target position, via \code{\link{idx_step_add}}. This is required
#' because heterogeneous steps (e.g. days mixed with months) cannot be
#' collapsed into a single scalar multiplier the way a uniform-unit
#' pattern can - so this constructor's calendars are O(distance from
#' anchor) to convert per position, rather than O(1).
#'
#' @param anchor A single reference point in calendar time (see
#' \code{\link{idx_calendar}}).
#' @param multi_step A \code{\link{multi_step_pattern}} object: the
#' repeating cycle of \code{idx_step}s.
#' @param anchor_pos,posixct,anchor_name As in \code{\link{idx_calendar}}.
#' @param pattern_start A single positive integer: the index into
#' \code{multi_step} at which \code{anchor_pos} falls (the slot for the
#' step *out of* \code{anchor_pos}). Defaults to \code{1}.
#'
#' @returns An object of class \code{idx_calendar}.
#'
#' @examples
#' # 3 days, 3 days, then 1 month, repeating.
#' cal <- idx_calendar_multi_step(
#'   anchor = as.Date("2024-01-01"),
#'   multi_step = multi_step_pattern(idx_step(days = 3), idx_step(days = 3),
#'                                    idx_step(months = 1)),
#'   posixct = TRUE
#' )
#' idx_to_date(cal, 1:5)
#'
#' @export
idx_calendar_multi_step <- function(anchor, multi_step, anchor_pos = 1L,
                                    pattern_start = 1L, posixct = FALSE,
                                    anchor_name = NULL) {
  if (!is_multi_step_pattern(multi_step)) {
    stop("multi_step must be a multi_step_pattern object; see multi_step_pattern().")
  }
  idx_calendar_validate_common(anchor_name, anchor_pos, pattern_start, posixct,
                               cycle_len = length(multi_step))
  structure(
    list(
      anchor = anchor,
      anchor_pos = as.integer(anchor_pos),
      amount = NA_real_,
      unit = NA_character_,
      step = NULL,
      multi_step = multi_step,
      pattern = NA_real_,
      pattern_start = as.integer(pattern_start),
      posixct = posixct,
      anchor_name = anchor_name
    ),
    class = "idx_calendar"
  )
}

#' @title Test whether an object is an \code{idx_calendar}
#' @param x Object to test.
#' @returns Logical.
#' @export
is_idx_calendar <- function(x) {
  inherits(x, "idx_calendar")
}

#' @title Print an \code{idx_calendar}
#' @param x An \code{idx_calendar} object.
#' @param ... Unused.
#' @export
print.idx_calendar <- function(x, ...) {
  cat("<idx_calendar>\n")
  if (!is.null(x$anchor_name)) cat("  anchor_name:   ", x$anchor_name, "\n")
  cat("  anchor:        ", format(x$anchor), "\n")
  cat("  anchor_pos:    ", x$anchor_pos, "\n")
  if (!is.null(x$multi_step)) {
    cat("  step:          multi_step_pattern (see below)\n")
  } else if (is.na(x$amount)) {
    cat("  step:          ")
    print(x$step)
  } else {
    cat("  step:          ", x$amount, x$unit, "\n")
  }
  if (!is.null(x$multi_step)) {
    print(x$multi_step)
  } else {
    cat("  pattern:       ", paste(x$pattern, collapse = " "), "\n")
  }
  cat("  pattern_start: ", x$pattern_start, "\n")
  cat("  posixct:       ", x$posixct, "\n")
  invisible(x)
}

#' @title Cumulative pattern-weighted step offset between two positions
#'
#' @description Computes, for one or more integer \code{idx_series}
#' positions, the number of base \code{amount}-sized steps that separate
#' each position from \code{cal$anchor_pos}, accounting for the repeating
#' \code{pattern} of step multipliers. This is a pure integer/numeric
#' computation - no calendar arithmetic and no units are applied here; the
#' result is in multiples of \code{cal$amount}, positive for positions
#' after the anchor and negative for positions before it. It is intended to
#' be consumed by \code{\link{idx_to_date}} (or directly by the plotting
#' layer's own POSIXct-aware conversion).
#'
#' @param cal An \code{idx_calendar} object.
#' @param pos An integer vector of \code{idx_series} positions.
#' @returns A numeric vector, the same length as \code{pos}, giving the
#' pattern-weighted step offset (in units of \code{cal$amount}) of each
#' position relative to \code{cal$anchor_pos}.
#' @export
idx_calendar_offset <- function(cal, pos) {
  stopifnot(is_idx_calendar(cal))
  pos <- as.integer(pos)
  pat <- cal$pattern
  plen <- length(pat)
  
  cum_from_start <- cumsum(pat)
  cycle_total <- cum_from_start[plen]
  
  weight_between <- function(a, b) {
    if (b >= a) {
      if (a == 1) cum_from_start[b] else cum_from_start[b] - cum_from_start[a - 1]
    } else {
      tail_w <- if (a == 1) 0 else cycle_total - cum_from_start[a - 1]
      tail_w + cum_from_start[b]
    }
  }
  
  offset <- numeric(length(pos))
  for (i in seq_along(pos)) {
    delta_pos <- pos[i] - cal$anchor_pos
    if (delta_pos == 0) {
      offset[i] <- 0
      next
    }
    fwd <- delta_pos > 0
    n_steps <- abs(delta_pos)
    n_full <- n_steps %/% plen
    r <- n_steps %% plen
    
    w <- n_full * cycle_total
    if (r > 0) {
      if (fwd) {
        start_slot <- cal$pattern_start
        end_slot <- ((start_slot - 1 + r - 1) %% plen) + 1
        w <- w + weight_between(start_slot, end_slot)
      } else {
        start_slot <- ((cal$pattern_start - 1 - r) %% plen) + 1
        end_slot <- ((cal$pattern_start - 1 - 1) %% plen) + 1
        w <- w + weight_between(start_slot, end_slot)
      }
    }
    offset[i] <- if (fwd) w else -w
  }
  offset
}

#' @title Translate positions to dates for a \code{multi_step_pattern}
#' calendar
#'
#' @description Internal helper for \code{\link{idx_to_date}}. Unlike the
#' numeric \code{pattern} mechanism, a \code{multi_step_pattern} cycles
#' through distinct \code{idx_step}s that generally cannot be collapsed
#' into a single scalar (e.g. mixing fixed-duration and calendar-relative
#' components), so each requested position is reached by walking one slot
#' at a time from \code{cal$anchor_pos}, applying \code{\link{idx_step_add}}
#' for the appropriate slot at each hop. This is O(distance from anchor)
#' per position, rather than the O(1) closed-form used for a numeric
#' \code{pattern}.
#'
#' @param cal An \code{idx_calendar} with a non-\code{NULL} \code{cal$multi_step}.
#' @param pos An integer vector of \code{idx_series} positions.
#' @returns A vector, the same length as \code{pos}, of calendar times.
#' @keywords internal
#' @noRd
idx_multi_step_offset_to_date <- function(cal, pos) {
  pos <- as.integer(pos)
  ms <- cal$multi_step
  plen <- length(ms)
  
  anchor0 <- if (isTRUE(cal$posixct)) as.POSIXct(cal$anchor) else cal$anchor
  
  walk_to <- function(target) {
    delta <- target - cal$anchor_pos
    if (delta == 0) return(anchor0)
    fwd <- delta > 0
    n_steps <- abs(delta)
    a <- anchor0
    slot <- cal$pattern_start
    for (i in seq_len(n_steps)) {
      if (fwd) {
        step_i <- ms[[slot]]
        a <- idx_step_add(a, step_i, 1)
        slot <- (slot %% plen) + 1
      } else {
        slot <- ((slot - 2) %% plen) + 1
        step_i <- ms[[slot]]
        a <- idx_step_add(a, step_i, -1)
      }
    }
    a
  }
  
  raw <- vapply(pos, function(p) unclass(walk_to(p)), FUN.VALUE = numeric(1))
  out <- raw
  attributes(out) <- attributes(anchor0)
  out
}

#' @title Translate a calendar date to an \code{idx_series} integer position,
#' for a \code{multi_step_pattern} calendar
#'
#' @description Internal helper; the inverse of
#' \code{\link{idx_multi_step_offset_to_date}} (in turn used by
#' \code{\link{idx_to_pos}} for \code{multi_step_pattern} calendars). Since a
#' \code{multi_step_pattern}'s slots generally cannot be collapsed into a
#' single scalar step size, there is no closed-form inverse either: this
#' walks one slot at a time from \code{cal$anchor_pos}, in whichever
#' direction \code{date} lies, comparing the resulting date to \code{date}
#' at each hop. Every \code{idx_step} component is constrained to be
#' non-negative (see \code{\link{idx_step}}) and at least one component of
#' each step is non-zero, so calendar time strictly increases with each
#' forward hop (strictly decreases with each backward hop); this
#' monotonicity is what makes the one-directional walk safe (no risk of
#' overshoot-and-return oscillation). This is O(distance from anchor), like
#' its forward counterpart.
#'
#' @param cal An \code{idx_calendar} with a non-\code{NULL} \code{cal$multi_step}.
#' @param date A single date/time, coerced to \code{cal$anchor}'s class.
#' @returns A single integer: the \code{idx_series} position corresponding
#' to \code{date}.
#' @keywords internal
#' @noRd
idx_multi_step_date_to_offset <- function(cal, date) {
  ms <- cal$multi_step
  plen <- length(ms)
  
  kind <- idx_anchor_kind(cal$anchor)
  
  date <- if (kind$is_yearqtr) {
    zoo::as.yearqtr(date)
  } else if (kind$is_yearmon) {
    zoo::as.yearmon(date)
  } else if (kind$is_posixct) {
    as.POSIXct(date, tz = format(cal$anchor, "%Z"))
  } else {
    as.Date(date)
  }

  as_num <- function(x) as.numeric(x)
  target_num <- as_num(date)
  anchor_num <- as_num(cal$anchor)
  
  if (isTRUE(all.equal(target_num, anchor_num))) {
    return(cal$anchor_pos)
  }
  fwd <- target_num > anchor_num
  
  a <- cal$anchor
  slot <- cal$pattern_start
  pos <- cal$anchor_pos
  # Generous but finite cap, so a date that never lands on a step boundary
  # (e.g. a typo, or a date on the wrong side of an irregular pattern)
  # fails with a clear error rather than looping indefinitely.
  max_hops <- 1e6L
  
  for (i in seq_len(max_hops)) {
    if (fwd) {
      step_i <- ms[[slot]]
      a_next <- idx_step_add(a, step_i, 1)
      a <- a_next
      pos <- pos + 1L
      slot <- (slot %% plen) + 1L
    } else {
      slot <- ((slot - 2) %% plen) + 1L
      step_i <- ms[[slot]]
      a_next <- idx_step_add(a, step_i, -1)
      a <- a_next
      pos <- pos - 1L
    }
    a_num <- as_num(a)
    if (isTRUE(all.equal(a_num, target_num))) {
      return(as.integer(pos))
    }
    if ((fwd && a_num > target_num) || (!fwd && a_num < target_num)) {
      stop("idx_to_pos: date does not fall on a whole step boundary of ",
           "cal's multi_step_pattern.")
    }
  }
  stop("idx_to_pos: date not reached within ", max_hops,
       " steps of cal$anchor_pos; check that date and cal$anchor use ",
       "compatible units/direction.")
}

#' @title Translate a calendar date to an \code{idx_series} integer position,
#' for a single compound \code{idx_step} calendar
#'
#' @description Internal helper; the inverse of \code{\link{idx_step_add}}
#' as used by \code{\link{idx_to_date}} for \code{\link{idx_calendar_step}}
#' calendars (a single compound \code{idx_step} repeated every position,
#' i.e. \code{cal$multi_step} is \code{NULL} but \code{cal$amount} is
#' \code{NA}). Since every \code{idx_step} field is non-negative and at
#' least one is non-zero, \code{idx_step_add(cal$anchor, cal$step, n)} is
#' strictly monotonic in \code{n}; this is what allows a direct binary
#' search on \code{n} rather than the slot-by-slot linear walk needed for
#' the general \code{multi_step_pattern} case (see
#' \code{\link{idx_multi_step_date_to_offset}}) - O(log distance) instead
#' of O(distance).
#'
#' @param cal An \code{idx_calendar} with \code{cal$multi_step} \code{NULL}
#' and \code{cal$amount} \code{NA} (i.e. built via
#' \code{\link{idx_calendar_step}}).
#' @param date A single date/time, coerced to \code{cal$anchor}'s class.
#' @returns A single integer: the \code{idx_series} position corresponding
#' to \code{date}.
#' @keywords internal
#' @noRd
idx_step_date_to_offset <- function(cal, date) {
  step <- cal$step
  kind <- idx_anchor_kind(cal$anchor)
  
  date <- if (kind$is_yearqtr) {
    zoo::as.yearqtr(date)
  } else if (kind$is_yearmon) {
    zoo::as.yearmon(date)
  } else if (kind$is_posixct) {
    as.POSIXct(date, tz = format(cal$anchor, "%Z"))
  } else {
    as.Date(date)
  }
  
  target_num <- as.numeric(date)
  at_n <- function(n) as.numeric(idx_step_add(cal$anchor, step, n))
  
  anchor_num <- as.numeric(cal$anchor)
  if (isTRUE(all.equal(target_num, anchor_num))) {
    return(cal$anchor_pos)
  }
  fwd <- target_num > anchor_num
  
  lo <- 0L; hi <- 1L
  repeat {
    probe <- if (fwd) hi else -hi
    v <- at_n(probe)
    reached_or_passed <- if (fwd) v >= target_num else v <= target_num
    if (reached_or_passed) break
    lo <- hi
    hi <- hi * 2L
    if (hi > 2e9) {
      stop("idx_to_pos: date is too far from cal$anchor to bracket via ",
           "cal$step (or does not fall on a step boundary in that direction).")
    }
  }
  lo_n <- if (fwd) lo else -lo
  hi_n <- if (fwd) hi else -hi
  
  while (abs(hi_n - lo_n) > 1) {
    mid_n <- lo_n + (hi_n - lo_n) %/% 2L
    v <- at_n(mid_n)
    on_target_side <- if (fwd) v < target_num else v > target_num
    if (on_target_side) lo_n <- mid_n else hi_n <- mid_n
  }
  for (n in c(lo_n, hi_n)) {
    if (isTRUE(all.equal(at_n(n), target_num))) {
      return(as.integer(cal$anchor_pos + n))
    }
  }
  stop("idx_to_pos: date does not fall on a whole step boundary of ",
       "cal$step relative to cal$anchor.")
}

#' @title Translate a date to a position for a numeric-pattern calendar
#'
#' @description Internal inverse of \code{\link{idx_to_date}} for calendars
#' built by \code{\link{idx_calendar}}. A numeric pattern makes the calendar
#' offset a piecewise sequence rather than a single amount that can be divided
#' in closed form, so this brackets and binary-searches the integer position.
#' Using \code{idx_to_date} as the forward map also verifies exact step
#' boundaries instead of silently rounding off-boundary dates.
#'
#' @param cal An \code{idx_calendar} built by \code{idx_calendar}.
#' @param date A single date/time.
#' @returns A single integer position.
#' @keywords internal
#' @noRd
idx_calendar_date_to_pos <- function(cal, date) {
  kind <- idx_anchor_kind(cal$anchor)
  date <- if (kind$is_yearqtr) {
    zoo::as.yearqtr(date)
  } else if (kind$is_yearmon) {
    zoo::as.yearmon(date)
  } else if (kind$is_posixct) {
    as.POSIXct(date, tz = format(cal$anchor, "%Z"))
  } else {
    as.Date(date)
  }
  
  calendar_month_granularity <- kind$is_date_posix &&
    cal$unit %in% c("months", "quarters", "years")
  comparable <- function(x) {
    if (calendar_month_granularity) {
      x <- as.Date(x)
      return(as.integer(format(x, "%Y")) * 12L +
               as.integer(format(x, "%m")) - 1L)
    }
    as.numeric(x)
  }
  target_num <- comparable(date)
  anchor_num <- comparable(cal$anchor)
  if (isTRUE(all.equal(target_num, anchor_num))) {
    return(as.integer(cal$anchor_pos))
  }
  fwd <- target_num > anchor_num
  at_delta <- function(delta) {
    comparable(idx_to_date(cal, cal$anchor_pos + delta))
  }
  
  lo <- 0
  hi <- 1
  repeat {
    probe <- if (fwd) hi else -hi
    value <- at_delta(probe)
    reached_or_passed <- if (fwd) value >= target_num else value <= target_num
    if (reached_or_passed) break
    lo <- hi
    if (hi > 1e9) {
      stop("idx_to_pos: date is too far from cal$anchor to bracket.")
    }
    hi <- hi * 2
  }
  lo_delta <- if (fwd) lo else -lo
  hi_delta <- if (fwd) hi else -hi
  
  while (abs(hi_delta - lo_delta) > 1) {
    mid_delta <- lo_delta + (hi_delta - lo_delta) %/% 2
    value <- at_delta(mid_delta)
    before_target <- if (fwd) value < target_num else value > target_num
    if (before_target) lo_delta <- mid_delta else hi_delta <- mid_delta
  }
  for (delta in c(lo_delta, hi_delta)) {
    if (isTRUE(all.equal(at_delta(delta), target_num))) {
      return(as.integer(cal$anchor_pos + delta))
    }
  }
  stop("idx_to_pos: date does not fall on a whole step boundary of ",
       "cal's amount/unit/pattern relative to cal$anchor.")
}

#' @title Translate \code{idx_series} positions to calendar time
#'
#' @description Converts one or more integer \code{idx_series} positions
#' into calendar time using an \code{idx_calendar}.
#'
#' When \code{cal$posixct} is \code{TRUE} and \code{cal$anchor} is a
#' \code{Date}/\code{POSIXct}, and \code{cal$unit} is one of the standard
#' calendar units - \code{"seconds"}, \code{"minutes"}, \code{"hours"},
#' \code{"days"}, \code{"weeks"}, \code{"months"}, \code{"quarters"},
#' \code{"years"} (singular or plural) - calendar-aware stepping is used via
#' \code{\link{seq.Date}}/\code{\link{seq.POSIXt}}'s \code{by} argument.
#' This is exact for every one of these units: seconds through weeks are
#' fixed-duration and so plain offset arithmetic would already agree with
#' it, but months/quarters/years are calendar-relative (variable length -
#' e.g. not every month has the same number of days) and genuinely require
#' this calendar-aware stepping to land on the correct date (respecting
#' month-end/DST rules) rather than a fixed number of days.
#'
#' \code{cal$anchor} may also be a \code{zoo} \code{yearqtr} or
#' \code{yearmon}, for series that are naturally indexed at quarterly or
#' monthly resolution rather than daily.
#'
#' Internally, every \code{idx_calendar} (whether built via the primary
#' \code{amount}/\code{unit} constructor or via
#' \code{\link{idx_calendar_step}}) carries a compound \code{cal$step}
#' (see \code{\link{idx_step}}), and calendar-aware stepping is delegated
#' uniformly to \code{\link{idx_step_add}}. This is what allows steps that
#' mix calendar-relative and fixed-duration components, e.g. "3 seconds
#' after each quarter" or "a month and 13 seconds" - a single
#' \code{amount}/\code{unit} pair cannot express these, but a compound
#' \code{idx_step} can.
#'
#' Otherwise (\code{posixct = FALSE}, or a non-Date/POSIXct/yearqtr/yearmon
#' anchor) the conversion falls back to plain numeric offset arithmetic: it
#' multiplies the pattern-weighted step offset (see
#' \code{\link{idx_calendar_offset}}) by \code{cal$amount} and adds the
#' result to \code{cal$anchor}. This supports arbitrary non-calendar units
#' (e.g. picoseconds) uniformly, but is not available for calendars built
#' via \code{\link{idx_calendar_step}} with \code{posixct = FALSE}, since
#' those have no single \code{amount} to multiply by (see
#' \code{\link{idx_calendar_step}}); such calendars require
#' \code{posixct = TRUE}.
#'
#' @param cal An \code{idx_calendar} object.
#' @param pos An integer vector of \code{idx_series} positions.
#' @returns A vector, the same length as \code{pos}, giving the calendar
#' time (or plain number, for non-calendar \code{cal$anchor}) for each
#' position.
#' @export
idx_to_date <- function(cal, pos) {
  stopifnot(is_idx_calendar(cal))
  
  is_calendar_anchor <- idx_anchor_kind(cal$anchor)$is_calendar
  
  if (!is.null(cal$multi_step)) {
    if (!isTRUE(cal$posixct) || !is_calendar_anchor) {
      stop("idx_to_date: a multi_step_pattern calendar requires posixct = TRUE ",
           "and a Date/POSIXct/yearqtr/yearmon anchor.")
    }
    return(idx_multi_step_offset_to_date(cal, pos))
  }
  
  offset <- idx_calendar_offset(cal, pos)
  
  if (isTRUE(cal$posixct) && is_calendar_anchor && !is.null(cal$step)) {
    if (!isTRUE(all.equal(offset, round(offset)))) {
      stop("idx_to_date: non-integer number of steps is not supported ",
           "for calendar-aware (posixct) offsets.")
    }
    idx_step_add(cal$anchor, cal$step, offset)
  } else {
    if (is.na(cal$amount)) {
      stop("idx_to_date: this idx_calendar was built via idx_calendar_step() ",
           "and has no single 'amount' to fall back to; posixct = TRUE is ",
           "required (and cal$anchor must be Date/POSIXct/yearqtr/yearmon).")
    }
    cal$anchor + offset * cal$amount
  }
}

#' @title Map an \code{idx_calendar} unit string to a \code{seq.Date}/
#' \code{seq.POSIXt} \code{by} unit
#'
#' @description Internal helper. Recognises the eight standard calendar
#' units that \code{POSIXct}/\code{Date} arithmetic can represent exactly -
#' seconds, minutes, hours, days, weeks, months, quarters, years (singular
#' or plural, case-insensitive) - and maps each to the singular form
#' expected by \code{seq.Date}/\code{seq.POSIXt}'s \code{by} argument.
#' Returns \code{NULL} for any other unit string, signalling that
#' calendar-aware stepping is not applicable and the caller should fall
#' back to plain numeric offset arithmetic.
#'
#' @param unit A single string (\code{cal$unit}).
#' @returns A single string suitable for \code{seq.Date(by = ...)}, or
#' \code{NULL} if \code{unit} is not one of the eight standard units.
#' @keywords internal
#' @noRd
idx_calendar_by_unit <- function(unit) {
  canonical <- idx_canonical_unit(unit)
  if (is.null(canonical)) NULL else canonical$singular
}

registerS3method("print", "idx_calendar", print.idx_calendar)
registerS3method("print", "idx_step", print.idx_step)
registerS3method("print", "multi_step_pattern", print.multi_step_pattern)