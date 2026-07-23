# Created by: Craig Thamotheram
# Created on: 15/02/2022
# Refactored: detach analysis from calendar time.

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' @title Index-anchored series
#'
#' @description \code{idx_series} is the core data structure used throughout
#' the analysis parts of this package in place of calendar-indexed objects
#' (such as \code{xts}). A series is stored as a plain numeric vector or
#' matrix of observations together with a single integer \code{start},
#' which is the integer position (relative to some external anchor point,
#' e.g. the first observation ever collected) of the first row of
#' \code{data}. All internal computation - estimation, filtering,
#' forecasting, reinitialisation - is done purely in terms of these integer
#' positions. This keeps the statistical engine of the package completely
#' free of calendar/time-of-day concerns.
#'
#' Calendar time is reintroduced only at the very edges of the package
#' (e.g. plotting, printing to the user) by a separate, cosmetic
#' translation layer that maps integer positions back to dates given an
#' anchor date. That translation layer is out of scope here.
#'
#' @param data A numeric vector or matrix of observations. For a matrix,
#' rows are observations (in position order) and columns are variables.
#' @param start A single positive integer giving the position of the first
#' observation in \code{data}. Defaults to \code{1}.
#'
#' @returns An object of class \code{idx_series}.
#'
#' @examples
#' x <- idx_series(cumsum(rpois(50, 5)) + 1, start = 1)
#' length(x)
#' idx_range(x)
#' x[10:20]
#'
#' @export
idx_series <- function(data, start = 1L) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (!is.numeric(data) && !is.matrix(data)) {
    stop("data must be a numeric vector or matrix.")
  }
  if (length(start) != 1 || !isTRUE(all.equal(start, as.integer(start)))) {
    stop("start must be a single integer.")
  }
  structure(
    list(data = data, start = as.integer(start)),
    class = "idx_series"
  )
}

# Register idx_series as an S4-compatible "old" S3 class so that it can be
# used as a field type in setRefClass()-based classes elsewhere in the
# package (mirrors setOldClass("KFS") for KFAS's KFS objects).
methods::setOldClass("idx_series")

#' @title Test whether an object is an \code{idx_series}
#' @param x Object to test.
#' @returns Logical.
#' @export
is_idx_series <- function(x) {
  inherits(x, "idx_series")
}

#' @title Coerce to \code{idx_series}
#'
#' @description Coerces plain vectors/matrices to \code{idx_series}, or
#' passes through an \code{idx_series} unchanged. Convenience helper so
#' that internal functions can accept either raw data or an already
#' constructed \code{idx_series}.
#'
#' @param x A numeric vector, matrix or \code{idx_series} object.
#' @param start Integer start position, only used if \code{x} is not
#' already an \code{idx_series}. Defaults to \code{1}.
#' @returns An \code{idx_series} object.
#' @export
as_idx_series <- function(x, start = 1L) {
  if (is_idx_series(x)) {
    return(x)
  }
  if (is.null(x)) {
    return(NULL)
  }
  idx_series(x, start = start)
}

#' @title Number of observations in an \code{idx_series}
#' @param x An \code{idx_series} object.
#' @param ... Unused.
#' @export
length.idx_series <- function(x, ...) {
  if (is.matrix(x$data)) NROW(x$data) else length(x$data)
}

#' @title Number of columns in an \code{idx_series}
#'
#' @description \code{base::NCOL()} is not an S3 generic, so it will not
#' dispatch on \code{idx_series} objects. Use this function instead when
#' you need the number of columns/variables in an \code{idx_series}.
#'
#' @param x An \code{idx_series} object.
#' @export
idx_ncol <- function(x) {
  stopifnot(is_idx_series(x))
  if (is.matrix(x$data)) ncol(x$data) else 1L
}

#' @title Number of columns in an \code{idx_series} (\code{NCOL} method)
#'
#' @description Provided for convenience when calling \code{NCOL()}
#' directly on an \code{idx_series}; note that since \code{base::NCOL} is
#' not an S3 generic, this method is only reached if \code{NCOL} is called
#' as \code{NCOL.idx_series()} directly, or after \code{NCOL} has been
#' made generic (e.g. via \code{Matrix} or similar). Prefer
#' \code{\link{idx_ncol}} for portable code within this package.
#'
#' @param x An \code{idx_series} object.
#' @export
NCOL.idx_series <- function(x) {
  idx_ncol(x)
}

#' @title Integer position range spanned by an \code{idx_series}
#'
#' @description Returns the first and last integer position covered by
#' \code{x}, i.e. \code{c(x$start, x$start + length(x) - 1)}.
#'
#' @param x An \code{idx_series} object.
#' @returns An integer vector of length 2: \code{c(first, last)}.
#' @export
idx_range <- function(x) {
  stopifnot(is_idx_series(x))
  n <- length(x)
  c(x$start, x$start + n - 1L)
}

#' @title Integer positions of an \code{idx_series}
#'
#' @description Returns the full vector of integer positions associated
#' with each observation in \code{x}. This is the index-based analogue of
#' \code{zoo::index()}.
#'
#' @param x An \code{idx_series} object.
#' @returns An integer vector.
#' @export
idx_positions <- function(x) {
  stopifnot(is_idx_series(x))
  n <- length(x)
  seq.int(x$start, length.out = n)
}

#' @title Subset an \code{idx_series} by position
#'
#' @description Subsets an \code{idx_series} using integer positions (in
#' the same units as \code{x$start}), not row numbers within \code{data}.
#' For example, if \code{x$start == 10}, then \code{x[10:12]} returns the
#' first three observations of \code{x}.
#'
#' @param x An \code{idx_series} object.
#' @param i Integer vector of positions to select.
#' @param j Optional column selector (for matrix-valued series).
#' @param ... Unused.
#' @export
`[.idx_series` <- function(x, i, j, ...) {
  pos <- idx_positions(x)
  rows <- match(i, pos)
  if (anyNA(rows)) {
    stop("Requested position(s) fall outside the range of this idx_series.")
  }
  if (is.matrix(x$data)) {
    newdata <- if (missing(j)) x$data[rows, , drop = FALSE] else x$data[rows, j, drop = FALSE]
  } else {
    newdata <- x$data[rows]
  }
  idx_series(newdata, start = min(i))
}

#' @title Print an \code{idx_series}
#' @param x An \code{idx_series} object.
#' @param ... Unused.
#' @export
print.idx_series <- function(x, ...) {
  cat("<idx_series> start =", x$start, ", n =", length(x), "\n")
  print(x$data)
  invisible(x)
}

#' @title Coerce an \code{idx_series} to a plain numeric vector or matrix
#'
#' @description \code{as.numeric()} is a primitive generic in base R and,
#' unlike ordinary closures, does not always reliably dispatch S3 methods
#' registered only under the \code{as.numeric.*} name (this depends on how
#' the call reaches the internal C-level coercion). To be safe,
#' \code{as.double.idx_series} is also registered, and both are
#' additionally registered explicitly via \code{registerS3method()} so
#' dispatch works regardless of how this file is sourced/loaded.
#'
#' @param x An \code{idx_series} object.
#' @param ... Unused.
#' @export
as.numeric.idx_series <- function(x, ...) {
  as.numeric(x$data)
}

#' @rdname as.numeric.idx_series
#' @export
as.double.idx_series <- function(x, ...) {
  as.numeric(x$data)
}

#' @title Explicit numeric extraction for \code{idx_series}
#'
#' @description Portable, always-dispatches accessor for the raw numeric
#' data underlying an \code{idx_series}. Prefer this over \code{as.numeric()}
#' in package-internal code, since it does not depend on S3 dispatch of a
#' primitive generic.
#'
#' @param x An \code{idx_series} object.
#' @returns A plain numeric vector (or matrix, if \code{x} wraps a matrix).
#' @export
idx_values <- function(x) {
  stopifnot(is_idx_series(x))
  x$data
}

#' @title Coerce an \code{idx_series} to a matrix
#' @param x An \code{idx_series} object.
#' @param ... Unused.
#' @export
as.matrix.idx_series <- function(x, ...) {
  as.matrix(x$data)
}

#' @title Combine (column-bind) \code{idx_series} objects
#'
#' @description Column-binds one or more \code{idx_series} objects that
#' share the same \code{start} and length, analogous to \code{cbind} for
#' matrices/xts objects.
#'
#' @param ... \code{idx_series} objects to combine.
#' @returns A single \code{idx_series} with combined columns.
#' @export
idx_cbind <- function(...) {
  args <- list(...)
  args <- args[!vapply(args, is.null, logical(1))]
  if (length(args) == 0) return(NULL)
  starts <- vapply(args, function(a) a$start, integer(1))
  lens <- vapply(args, length, integer(1))
  if (length(unique(starts)) > 1 || length(unique(lens)) > 1) {
    stop("All idx_series objects must share the same start and length to be combined.")
  }
  mats <- lapply(args, function(a) as.matrix(a$data))
  idx_series(do.call(cbind, mats), start = starts[1])
}

#' @title Take a lagged difference of an \code{idx_series}
#'
#' @description Index-based analogue of \code{diff()} for calendar series:
#' returns \eqn{x_t - x_{t-lag}} for each position \eqn{t}, dropping the
#' first \code{lag} positions (which have no valid lag within the series).
#'
#' @param x An \code{idx_series} object.
#' @param lag Integer number of positions to lag by. Default \code{1}.
#' @returns An \code{idx_series}, shorter than \code{x} by \code{lag}.
#' @export
idx_diff <- function(x, lag = 1L) {
  stopifnot(is_idx_series(x))
  n <- length(x)
  if (lag >= n) stop("lag must be smaller than the length of x.")
  if (is.matrix(x$data)) {
    newdata <- x$data[(lag + 1):n, , drop = FALSE] - x$data[1:(n - lag), , drop = FALSE]
  } else {
    newdata <- x$data[(lag + 1):n] - x$data[1:(n - lag)]
  }
  idx_series(newdata, start = x$start + lag)
}

#' @title Lag an \code{idx_series}
#'
#' @description Index-based analogue of \code{stats::lag()}: shifts the
#' series forward by \code{k} positions. The returned series still covers
#' \code{length(x)} observations, but its \code{start} is shifted by
#' \code{k}, i.e. position \eqn{t} of the result holds the value that
#' was at position \eqn{t-k} in \code{x}.
#'
#' @param x An \code{idx_series} object.
#' @param k Integer number of positions to lag by. Default \code{1}.
#' @returns An \code{idx_series} with the same length as \code{x}.
#' @export
idx_lag <- function(x, k = 1L) {
  stopifnot(is_idx_series(x))
  idx_series(x$data, start = x$start + k)
}

#' @title Bind two \code{idx_series} sequentially
#'
#' @description Concatenates \code{x} then \code{y} in position order.
#' \code{y} must start exactly one position after \code{x} ends.
#'
#' @param x,y \code{idx_series} objects to concatenate.
#' @returns A single, combined \code{idx_series}.
#' @export
idx_rbind <- function(x, y) {
  stopifnot(is_idx_series(x), is_idx_series(y))
  if (y$start != idx_range(x)[2] + 1L) {
    stop("y must start immediately after x ends.")
  }
  if (is.matrix(x$data) || is.matrix(y$data)) {
    newdata <- rbind(as.matrix(x$data), as.matrix(y$data))
  } else {
    newdata <- c(x$data, y$data)
  }
  idx_series(newdata, start = x$start)
}

## ----------------------------------------------------------------------
## Explicit S3 method registration.
##
## When this package is loaded normally (library(tsgc)), methods declared
## with @export/@S3method in NAMESPACE are registered automatically. When
## this file is merely source()'d (e.g. in ad hoc scripts or this file's
## own test harness), that registration can be skipped for some primitive
## generics (notably as.numeric/as.double), causing dispatch to silently
## fall through to the default method instead of *.idx_series. Registering
## explicitly here makes behaviour identical in both loading paths.
## ----------------------------------------------------------------------
registerS3method("length", "idx_series", length.idx_series)
registerS3method("[", "idx_series", `[.idx_series`)
registerS3method("print", "idx_series", print.idx_series)
registerS3method("as.numeric", "idx_series", as.numeric.idx_series)
registerS3method("as.double", "idx_series", as.double.idx_series)
registerS3method("as.matrix", "idx_series", as.matrix.idx_series)
# Note: base::NCOL() is a plain function (it does not call UseMethod()),
# so it can never be made to dispatch on idx_series via S3 registration.
# Use idx_ncol() instead of NCOL() for idx_series objects throughout this
# package.