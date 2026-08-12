# tsgc 2.0.0

## Breaking changes

* Internal computation no longer relies on `xts`/`zoo`. Series used during
  estimation and filtering are now represented with a new, dependency-free
  `idx_series` class, with calendar handling provided by a companion
  `idx_calendar` class. Bundled datasets are unaffected and remain `xts`
  objects; `xts_to_idx()` is provided to convert user-supplied `xts`/`zoo`
  data to `idx_series` where needed. Code that relied on the internal
  (non-exported) representation of series may need to be updated.

* `SSModelDynGompertzReinit` has been deprecated. Its functionality is now
  provided directly by `SSModelDynamicGompertz` via the `reinit.idx`,
  `original.results`, and `use.presample.info` fields, so a single class
  now covers both the standard and reinitialised model. `SSModelBase` has
  also been deprecated, with its logic folded into `SSModelDynamicGompertz`.

* In `FilterResults`, the `confidence_level` argument to `predict_all()`
  and `get_gy_ci()` has been renamed to `confidence.level`, for
  consistency with naming elsewhere in the package. Calls that pass this
  argument by name will need to be updated. The default value of
  `sea.on` in `predict_all()` has also changed, from `FALSE` to `TRUE`.

* `forecast_peak()`, `forecast.peak()`, and `plot_new_cases()` have been
  deprecated.

## New features

* Added a leading indicator model, via the new `SSModelLeadingIndicator`
  and `FilterResultsLI` reference classes, for modelling a target series
  using a leading series (see Harvey, A. (2021). "Time Series Modelling
  of Epidemics: Leading Indicators, Control Groups and Policy
  Assessment." *National Institute Economic Review*, 257, 83-100).
  `FilterResultsLI` provides `predict_level()`, `predict_all()`,
  `get_growth_y()`, `get_gy_ci()`, `print()`, `summary()`, and `mapes()`
  methods, paralleling `FilterResults` for the standard model.

* Added five new example datasets: `ukitaly`, `gauteng_weather_2021`,
  `england_weather_2021`, `nintendo_sales`, and `etrading_apps`, including
  weather data to support leading indicator examples.

* Added `idx_series()` and supporting methods (`idx_cbind()`, `idx_rbind()`,
  `idx_diff()`, `idx_lag()`, `idx_range()`, `idx_values()`,
  `idx_positions()`, and standard methods including `length()`, `head()`,
  `tail()`, `as.matrix()`, `as.numeric()`, `as.double()`, and `print()`)
  as a dependency-free representation of time series for internal
  computation.

* Added `idx_calendar()` and supporting helpers (`idx_step()`,
  `idx_to_date()`, `idx_calendar_offset()`, `idx_calendar_step()`,
  `idx_calendar_multi_step()`, and related utilities) providing a
  consistent translation layer between integer positions and calendar
  dates for daily, weekly, monthly, quarterly, and annual data.

* Added `get_timeframe()`, `idx_to_pos()`, `idx_offset_to_pos()`, and
  `idx_detect_calendar_pattern()` to support the new indexing system.

* Added `check_variance_boundary()` and `print_model_diagnostics()` to
  surface variance boundary conditions and model diagnostics more clearly.

* Added `xts_to_idx()` to convert `xts`/`zoo` objects to the internal
  `idx_series` representation, for users supplying their own `xts` data.

* Added `cross_val()` for cross-validating one or more models over a
  series of rolling estimation windows, and `mapes()` for computing mean
  absolute percentage errors from forecast results.

* Added `df2ldl_lead()` for handling two-column data for leading 
  indicator models.

* Added `print()` and `summary()` S3 methods for `SSModelDynamicGompertz`,
  `SSModelLeadingIndicator`, `FilterResults`, and `FilterResultsLI`.

* Added new plotting functions `plot_compare_forecast()` (compare
  forecasts across models), `plot_log_forecast()` (forecasts on the log
  scale), and `plot_r0()` (plot estimated reproduction numbers from
  `estimate_r0()`). Also added plotting functionality for models 
  (`plot()` takes in a model object)

* Added a new set of accessor functions (`accessorFns.R`) for extracting
  components from fitted model and `KFS` objects, including `output()`,
  `modelKFS()`, `seasonalComp()`, `att()`, `Ptt()`, `get_V()`, `gety()`,
  `alphahat()`, and a standalone `estimate()` wrapper around a model's
  `estimate()` method.

## Improvements

* Substantially expanded plotting internals with helper functions
  (`idx_series_df()`, `idx_resolve_axis()`, `idx_axis_opts()`,
  `idx_x_scale()`, `idx_x_lab()`, `idx_add_info_box()`, and related
  helpers) to give more consistent, better-labelled axes across all plot
  types, including support for non-daily frequencies.

* Reduced the package's dependency footprint by removing reliance on
  `xts`/`zoo` from internal computation.

## Bug fixes

* Various small fixes to filtering, forecasting, and plotting uncovered
  during the internal rewrite.

# tsgc 0.0 (2022-02-08)

* Added a `NEWS.md` file to track changes to the package.
