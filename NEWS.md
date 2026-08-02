# tsgc 2.0.0

## Breaking changes

* The internal computational engine no longer relies on `xts`/`zoo` for
  representing time series during estimation and filtering. Series are now
  represented internally with a new, lightweight `idx_series` class, and
  calendar handling is provided by a companion `idx_calendar` class. This is
  purely an internal architectural change for most users, but code that
  relied on the internal (non-exported) representation of series may need
  to be updated.
* The `plot()` method has been removed from the `SSModelDynamicGompertz`
  and `SSModelLeadingIndicator` R6 classes. Plotting is now handled
  exclusively via the S3 methods `plot.SSModelDynamicGompertz()` and
  `plot.SSModelLeadingIndicator()`, which continue to work with the usual
  `plot(model, ...)` syntax.
* The `timetk` package is no longer a dependency.
* `estimate_r0()` now returns a data frame, including dates from
  `res$calendar` when available. The redundant `estimate_r0_df()` wrapper
  has been removed.

## New features

* Added `idx_series()` and supporting methods (`idx_cbind()`, `idx_rbind()`,
  `idx_diff()`, `idx_lag()`, `idx_range()`, `idx_values()`,
  `idx_positions()`, and standard methods such as `length()`, `head()`,
  `tail()`, `as.matrix()`, `as.numeric()`, `as.double()`, and `print()`) as
  a dependency-free replacement for calendar-indexed objects in internal
  computation.
* Added `idx_calendar()` and supporting helpers (`idx_step()`, `idx_to_date()`,
  `idx_calendar_offset()`, `idx_calendar_step()`, `idx_calendar_multi_step()`,
  and related utilities) providing a single, consistent translation layer
  between integer positions and calendar dates for daily, weekly, monthly,
  quarterly, and annual data.
* Added `get_time_resolution()`-style helpers `get_timeframe()`,
  `idx_to_pos()`, `idx_offset_to_pos()`, and
  `idx_detect_calendar_pattern()` to `utils.R` to support the new indexing
  system.
* Added `check_variance_boundary()` and `print_model_diagnostics()` to
  surface variance boundary conditions and model diagnostics more clearly
  to users.
* Added `xts_to_idx()` to convert `xts` objects into the internal
  `idx_series` representation, easing interoperability for users supplying
  `xts` data.

## Improvements

* Substantially expanded plotting internals (`plotting.R`) with helper
  functions (`idx_series_df()`, `idx_resolve_axis()`, `idx_axis_opts()`,
  `idx_x_scale()`, `idx_x_lab()`, `idx_add_info_box()`, and related
  helpers) to give more consistent, better-labelled axes across all plot
  types, including support for non-daily frequencies.
* Reduced the package's dependency footprint by removing reliance on
  `xts`/`zoo` from internal computation and dropping the unused `timetk`
  dependency.
* Expanded and reorganised the test suite, including a dedicated
  boundary-conditions test file (`test-boundary-conditions.R`) and new
  tests for the `idx_series`/`idx_calendar` system.
* Updated the vignette and replication script to reflect the internal
  changes above; user-facing behaviour and function signatures for
  `SSModelDynamicGompertz`, `SSModelLeadingIndicator`, `cross_val()`,
  and the plotting functions are unchanged.

## Bug fixes

* Various small fixes to filtering, forecasting, and plotting uncovered by
  the expanded test suite (see GitHub commit history for details).

# tsgc 0.0 (2022-02-08)

* Added a `NEWS.md` file to track changes to the package.
