## New release

This is a release. It is a major version update (from 0.0 to 2.0.0) that
includes breaking changes. This new release also includes a change of package
maintainer, from Craig Thamotheram <craig_thamotheram@hotmail.com> to 
Michael Ashby <mwa22@cam.ac.uk>. Craig Thamotheram has separately emailed
<CRAN-submissions@R-project.org> to confirm this change. Please note that
Craig's email address is no longer <cpt@tacindex.com> as he has since left
the company.

In this version we have:

* Rewritten the internal computational engine to use a new lightweight,
  dependency-free `idx_series`/`idx_calendar` indexing abstraction in place
  of `xts`/`zoo` for internal state-space calculations. `xts` and `zoo` are retained
  only at the user-facing edges (bundled datasets, `xts_to_idx()`, and
  plotting) where calendar-indexed objects are convenient.
* Added a leading indicator model (`SSModelLeadingIndicator` and
  `FilterResultsLI`), which models a target series using a related series
  that moves ahead of it, together with five new example datasets and
  supporting plotting functions.
* Consolidated `SSModelDynGompertzReinit` and `SSModelBase` into
  `SSModelDynamicGompertz`, so a single reference class now covers both
  the standard and reinitialised model.
* Moved the `plot()` method off the `SSModelDynamicGompertz` and
  `SSModelLeadingIndicator` reference classes and consolidated all
  plotting behaviour behind the corresponding S3 `plot.*` generics,
  avoiding duplicate plotting code paths.
* Added `kableExtra` to `Imports`, reflecting its use in the new
  `print_estimation_results()` output methods.
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