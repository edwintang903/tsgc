#' @title Calendar translation layer for \code{idx_series}
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
  if (!is.null(anchor_name) && (length(anchor_name) != 1 || !is.character(anchor_name))) {
    stop("anchor_name must be NULL or a single string.")
  }
  if (length(anchor_pos) != 1 || !isTRUE(all.equal(anchor_pos, as.integer(anchor_pos)))) {
    stop("anchor_pos must be a single integer.")
  }
  if (length(amount) != 1 || !is.numeric(amount) || !is.finite(amount) || amount <= 0) {
    stop("amount must be a single positive, finite number.")
  }
  if (length(unit) != 1 || !is.character(unit)) {
    stop("unit must be a single string.")
  }
  if (!is.numeric(pattern) || length(pattern) < 1 || any(!is.finite(pattern)) || any(pattern <= 0)) {
    stop("pattern must be a numeric vector of one or more positive, finite values.")
  }
  if (length(pattern_start) != 1 || !isTRUE(all.equal(pattern_start, as.integer(pattern_start))) ||
      pattern_start < 1 || pattern_start > length(pattern)) {
    stop("pattern_start must be a single integer between 1 and length(pattern).")
  }
  if (length(posixct) != 1 || !is.logical(posixct) || is.na(posixct)) {
    stop("posixct must be a single logical (TRUE/FALSE).")
  }
  structure(
    list(
      anchor = anchor,
      anchor_pos = as.integer(anchor_pos),
      amount = amount,
      unit = unit,
      pattern = as.numeric(pattern),
      pattern_start = as.integer(pattern_start),
      posixct = posixct,
      anchor_name = anchor_name
    ),
    class = "idx_calendar"
  )
}

methods::setOldClass("idx_calendar")

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
  cat("  step:          ", x$amount, x$unit, "\n")
  cat("  pattern:       ", paste(x$pattern, collapse = " "), "\n")
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
  
  # Cumulative pattern weight ending at each cycle slot.
  cum_from_start <- cumsum(pat)
  cycle_total <- cum_from_start[plen]
  
  # Pattern-weighted steps from cycle slot a to cycle slot b, moving
  # forward and wrapping around the pattern if necessary.
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
#' Otherwise (\code{posixct = FALSE}, a non-Date/POSIXct anchor, or a unit
#' outside the standard eight) the conversion falls back to plain numeric
#' offset arithmetic: it multiplies the pattern-weighted step offset (see
#' \code{\link{idx_calendar_offset}}) by \code{cal$amount} and adds the
#' result to \code{cal$anchor}. This supports arbitrary non-calendar units
#' (e.g. picoseconds) uniformly.
#'
#' @param cal An \code{idx_calendar} object.
#' @param pos An integer vector of \code{idx_series} positions.
#' @returns A vector, the same length as \code{pos}, giving the calendar
#' time (or plain number, for non-calendar \code{cal$anchor}) for each
#' position.
#' @export
idx_to_date <- function(cal, pos) {
  stopifnot(is_idx_calendar(cal))
  offset <- idx_calendar_offset(cal, pos)
  
  is_calendar_anchor <- inherits(cal$anchor, "Date") || inherits(cal$anchor, "POSIXct")
  by_unit <- idx_calendar_by_unit(cal$unit)
  
  if (isTRUE(cal$posixct) && is_calendar_anchor && !is.null(by_unit)) {
    step <- offset * cal$amount
    if (!isTRUE(all.equal(step, round(step)))) {
      stop("idx_to_date: non-integer number of ", cal$unit,
           " steps is not supported for calendar-aware (posixct) offsets.")
    }
    step <- as.integer(round(step))
    vapply(step, function(s) {
      seq(cal$anchor, by = paste(s, by_unit), length.out = 2)[2]
    }, FUN.VALUE = cal$anchor[1])
  } else {
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
  u <- tolower(unit)
  switch(u,
         second = , seconds = "sec",
         minute = , minutes = "min",
         hour = , hours = "hour",
         day = , days = "day",
         week = , weeks = "week",
         month = , months = "month",
         quarter = , quarters = "quarter",
         year = , years = "year",
         NULL
  )
}

registerS3method("print", "idx_calendar", print.idx_calendar)