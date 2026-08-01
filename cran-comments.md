## Resubmission

This is a resubmission. In this version I have:

* Rewritten the internal computational engine to use a new lightweight,
  dependency-free `idx_series`/`idx_calendar` indexing abstraction in place
  of `xts`/`zoo` for internal state-space calculations, substantially
  reducing the package's dependency footprint. `xts` and `zoo` are retained
  only at the user-facing edges (import/export and plotting) where
  calendar-indexed objects are convenient.
* Removed the unused `timetk` dependency from `Imports`.
* Moved the `plot()` method off the `SSModelDynamicGompertz` and
  `SSModelLeadingIndicator` R6 classes and consolidated all plotting
  behaviour behind the corresponding S3 `plot.*` generics, avoiding
  duplicate plotting code paths.
* Expanded unit test coverage, including new boundary-condition tests.
* Updated the vignette and replication script to reflect the above changes.

## Test environments

* local: macOS, R 4.4.x
* win-builder: R-devel, R-release, R-oldrelease
* R-hub: Windows Server 2022 (R-devel), Ubuntu Linux (R-release), Fedora
  Linux (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.
