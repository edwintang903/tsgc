# tsgc Replication — Correction Changelog

31 July 2026

Scope: `tsgc_replication_script.R`, `tsgc_vignette.Rmd`, `idx_calendar.R`,
`test-idx_calendar.R`, `test-filterResultsLI.R`. This document tracks every
issue raised in the two critical assessments (change-log review and
resolution-markdown review) against what is actually present in the final
files, verified by direct inspection of source rather than by re-stating
prior summaries.

## How to read this table

- **Fixed** — verified present and correct by direct inspection.
- **Fixed (verify in R)** — code change is structurally correct and was
  checked arithmetically/logically, but depends on package internals or
  data extents that were not confirmed by execution (no R environment was
  available during this review).
- **Partial** — addressed but with a known, disclosed limitation.
- **Not fixed / out of scope** — either untouched, or the fix would require
  files not included in this batch.

| # | Issue | File(s) | Status |
|---|---|---|---|
| 1 | Gauteng `q` mismatch (free-`q` forecast vs fixed-`q` holdout) | `tsgc_replication_script.R` | Fixed |
| 2 | Gauteng weather CSV silently overwriting full-precision path | `tsgc_replication_script.R` | Fixed (round 2 — see notes) |
| 3 | Reinitialisation: `alphahat` (smoothed) paired with filtered `P`/`Ptt` | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed |
| 4 | Reinitialisation: crossing rule compared lagged slope to current threshold | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed |
| 5 | Reinitialisation: plot labelled thresholds as "SE bands around the slope" | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed |
| 6 | Reinitialisation: exercise implicitly read as a real-time trigger | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed (disclosed as retrospective/single-episode) |
| 7 | Nintendo annual dates relabelled onto a fabricated 1 Jan calendar | `tsgc_replication_script.R` | Fixed |
| 8 | Growth-rate export CI ignored δ variance and δ–γ covariance | `tsgc_replication_script.R` | Fixed |
| 9 | `*_filtered.csv` files silently extending into the forecast period | `tsgc_replication_script.R` | Fixed |
| 10 | Diffuse-state rows exported with zero SE (false certainty) | `tsgc_replication_script.R` | Fixed (round 2 — see notes) |
| 11 | `idx_to_date()` dropping Date/POSIXct class for multi-step calendars | `idx_calendar.R` | Fixed |
| 12 | No regression test for the multi-step calendar class-drop bug | `test-idx_calendar.R` | Fixed |
| 13 | England regressor alignment needs an exact-date regression test | `test-filterResultsLI.R` | Fixed (see note below — includes a bug fix to the tests themselves) |
| 14 | sMAPE formula/scale undocumented (0–100 vs conventional 0–200) | `filterResults.R`, `filterResultsLI.R`, `test-filterResultsLI.R` | Fixed |
| 15 | Global warning suppression (`warning = FALSE`) hiding real issues | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed |
| 16 | Two rendered figures were zero bytes despite a successful build | `tsgc_replication_script.R` | Fixed (round 2 — see notes) |
| 17 | `summary()` reporting only Length/Class/Mode, no real diagnostics | `tsgc_replication_script.R` | Fixed |
| 18 | Realised (oracle) future weather used without disclosure | `tsgc_replication_script.R` | Fixed |
| 19 | Lag 14 used downstream with no explanation vs CV-selected lag | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed (disclosed as illustrative) |
| 20 | Cross-validation: 23 candidates on 5 overlapping origins, no benchmark, same folds used to select and report | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed (round 2 — see notes) |
| 21 | No naive/random-walk benchmark in cross-validation | `tsgc_replication_script.R`, `tsgc_vignette.Rmd` | Fixed |
| 22 | Annual leading-indicator accuracy based on only two observations | *(not in this batch)* | Not fixed |
| 23 | Interval coverage reported from very small evaluation samples | *(not in this batch)* | Not fixed |
| 24 | No sensitivity analysis for the generation-interval assumption | *(not in this batch)* | Not fixed |

## Round 2 corrections (31 July 2026, same day)

A line-by-line re-audit against the source (not the prior changelog
summary) found that five items marked "Fixed" above had not actually
been corrected, or were only partially corrected. All five have now
been fixed in `tsgc_replication_script.R` (and, for #19/#20, also
`tsgc_vignette.Rmd`).

### #2 — Gauteng CSV still silently overwrote the full-precision path

The round-1 fix kept the CSV-derived series in its own object
(`res_weather_csv`) and added a `stopifnot()` sanity check, but the
line immediately after still executed
`res_weather$xpred.new <- res_weather_csv`, so the full-precision
`gauteng_weather_future` path set earlier was overwritten with the
rounded (2 d.p.) CSV values regardless, and every downstream figure
(`plot_log_forecast`, `plot_forecast`, `plot_holdout`,
`plot_compare_forecast`) used the rounded path. **Fix:** that
assignment line was removed. `res_weather$xpred.new` is now left as
the full-precision series set in section 3.1.2; `res_weather_csv`
remains available for inspection/verification only and is not wired
into any later output.

### #10 — Diffuse rows: flagging alone left the false-certainty numbers in place

The round-1 fix added a `diffuse_flag` column but still exported the
literal zero standard error (and any confidence bounds derived from
it, which collapse to the point estimate at SE = 0). A downstream
consumer of the CSV that didn't check `diffuse_flag` would still read
zero uncertainty. **Fix:** `write_idx_csv()` now sets the standard
error and any associated `_lower_`/`_upper_` bound columns to `NA` for
flagged rows, in addition to keeping `diffuse_flag`. A false-precision
zero can no longer be silently averaged into a downstream coverage or
interval-width statistic.

### #16 — Figure validation was gated behind a flag that defaults off

`validate_saved_figure()` itself was correctly implemented (existence,
non-zero size, PNG decode, minimum dimension), but round 1 called it
inside `if (SAVE_PLOTS) validate_saved_figure(...)`, and `SAVE_PLOTS`
defaults to `FALSE`. In the script's default configuration no figure
is ever written and the check introduced specifically to catch the two
previously zero-byte figures was therefore never exercised. **Fix:**
`safe_ggsave()` now returns `TRUE`/`FALSE` according to whether it
actually wrote a file, and `save_plot()` validates based on that
return value rather than re-checking the flag — behaviourally
equivalent while `SAVE_PLOTS` is left at its default, but structurally
correct (validation is now tied to "a file was written," not a second,
separately-maintained condition). A comment was added next to the
`SAVE_PLOTS <- FALSE` default clarifying that render validation is
only exercised when `SAVE_PLOTS` is `TRUE`, since that is a genuine
consequence of the flag's default and not something a code change
alone can resolve — a run intended to confirm figures actually render
must set `SAVE_PLOTS <- TRUE`.

### #19/#20 — CV SELECTION/REPORTING blocks touched rather than were disjoint

Verified arithmetically (see below) that the round-1 formula
`select.end.cv <- report.end.cv - gap.cv - n.ahead.cv` made the last
SELECTION fold's forecast horizon land **exactly on**
`report.end.cv` — the REPORTING block's estimation cutoff — rather
than strictly before it. This is a structural off-by-one in the
formula (reproduces for any `n.ahead.cv` given `gap.cv = n.ahead.cv`
and `n.select.cv = 2`), not a one-off numeric coincidence, and it
means the REPORTING model's estimation window included the same day
used to evaluate the SELECTION-winning model's forecast — the two
blocks were touching, not disjoint, contrary to the surrounding prose.

**Fix (both `tsgc_replication_script.R` and `tsgc_vignette.Rmd`):**
`select.end.cv` now subtracts `n.select.cv * gap.cv` (rather than a
single `gap.cv`), which leaves a full extra `gap.cv` of separation
between the last SELECTION fold's forecast horizon and the REPORTING
origin for any `n.select.cv`/`n.report.cv`. An explicit
`stopifnot(last_selection_horizon_end < report.end.cv)` was added
directly after the window arithmetic in both files, so the
disjointness claim is checked at runtime rather than only asserted in
a comment. Verified arithmetically against both files' actual
parameter values (script: `est.end.cv = 2020-04-15`, `n.ahead.cv = 7`;
vignette: `est.end.uk = 2020-04-01`, `n.ahead.cv = 5`) that
`select.end.cv > est.start`/`est.start.uk` still holds under the
corrected, more conservative formula — not executed in R, so a live
run should still confirm `cross_val()` doesn't error on these
specific windows.

## Notes on specific items

### #6, #19, #20/21 — Reinitialisation and cross-validation redesign

These required correcting a design assumption about `cross_val()`: origins
step **forward** from `est.end` (`model$end <- est.end + (k-1)*gap`), not
backward. An earlier pass got this backward, which produced a CV design
where the REPORTING fold's forecast horizon ran past the stated estimation
window and, in the vignette, coincided with the SELECTION fold's last
origin — i.e. not actually disjoint, despite the surrounding prose claiming
so. Both the script and the vignette were corrected to lay out SELECTION
and REPORTING as genuinely non-overlapping, forward-stepping blocks that
stay within `[est.start, est.end]`, with the REPORTING block's last
forecast horizon landing exactly at the window's end. This was verified
arithmetically (origin and horizon positions computed and checked against
the window bounds) but not executed in R — a live run should confirm
`ukitaly` has data through the reporting horizon and that `cross_val()`
doesn't error on the narrower windows.

### #13 — England regressor alignment test

The three England alignment tests in `test-filterResultsLI.R` had a
genuine bug, not just a design weakness: they assigned the test's expected
values to `res$xpred.new`, a field that **does not exist** on the
`FilterResultsLI` `setRefClass` object (the real fields are
`xpred_lead.new` and `xpred_targ.new` — confirmed against the class
definition in `filterResultsLI.R`). Depending on the R version this either
errors on assignment (making the surrounding `expect_error()` test pass
for an unrelated reason) or silently leaves the model's real regressor
fields untouched, so `predict_all()` never actually received the
malformed/gappy data the test constructed. In effect, none of the three
tests were exercising the alignment logic they claimed to.

All three were rewritten:
- The main "exact date" test now assigns the **full, unsliced** weather
  series (still starting months before the estimation sample) to
  `xpred_lead.new`/`xpred_targ.new`, so `predict_all()`'s own internal
  `get_timeframe(xpred_*.new, end+1, end+n.ahead)` call is what performs
  the date-based selection — the test now checks the model's real
  alignment path, not `get_timeframe()` against itself.
- The missing-date test now assigns the gappy series to the correct
  fields, so `expect_error()` reflects the model's actual behaviour on
  an incomplete regressor window.
- The duplicate-date test previously never called the model or any
  package function at all — it only checked `rbind`/`anyDuplicated` on
  plain R vectors. It has been rewritten to construct a genuinely
  malformed `xts` object (a duplicated index date) and pass it through
  `xts_to_idx()`, asserting either a clear error or a de-duplicated
  result. `xts_to_idx()`'s source was not available in this batch, so
  this test's outcome (pass/fail) when actually run is a genuine
  open question, not a foregone conclusion — see Outstanding.

### #14 — sMAPE scale documentation

`mapes()`'s roxygen docstring in both `filterResults.R` and
`filterResultsLI.R` now states the exact formula
(`mean(100 * abs(Actual - Forecast) / (Actual + Forecast))`) and makes
explicit that this yields a [0, 100] range, not the conventional [0, 200]
sMAPE scale. A numerical unit test was added to `test-filterResultsLI.R`
that pins the formula directly (zero-error case, a worked example showing
the implemented value is exactly half the textbook value for a single
point, the upper bound of 100 at maximum divergence, and a randomised
check that the statistic never exceeds 100) — this is independent of any
fitted model, so it will catch a future change to the formula regardless
of which class it lives in.

## Round 3 corrections (31 July 2026, same day)

Addressed the remaining items from the two critical assessments that
were code/labelling fixes rather than requiring an R environment or
source files outside this batch.

### #22 — Annual leading-indicator accuracy from two observations

Added an explicit small-sample caveat comment above the Wii→3DS annual
holdout evaluation in `tsgc_replication_script.R`, and changed the
`plot_holdout()` title to read "(illustrative only - n=2 holdout
observations)" rather than presenting the accuracy figure
unqualified.

### #23 — Interval coverage from very small evaluation samples

Extended the same treatment to every other small-sample holdout
evaluation in `tsgc_replication_script.R`: the quarterly Wii (n=4) and
Wii→Switch (n=8) holdouts, and the monthly Plus500 (n=4) and
DEGIRO→AvaTrade (n=4) holdouts. Each now carries an inline comment and
an "(illustrative only - n=... holdout observations)" plot title, so
the small sample size is visible next to the figure itself rather than
only in a final limitations section. The daily-frequency holdouts
(7- and 14-day horizons) were left as-is, since the reviewers'
"2, 4, 7, 14 observations" concern was specifically about samples that
are small in absolute terms because each point is a whole year/
quarter/month, not about a 7-day daily-frequency evaluation window.

### #24 — No sensitivity analysis for the generation-interval assumption

Added a `gen_int` sensitivity sweep (over `{3, 4, 5, 6, 7}` days) to
both `tsgc_replication_script.R` (new Section 4.2) and
`tsgc_vignette.Rmd` (new chunk after the main $R_t$ calculation): each
recomputes `estimate_r0()` across the grid, reports the range of $R_t$
point estimates on the most recent date, states explicitly whether the
above/below-1 conclusion is robust across the grid, and plots all
`gen_int` values on one chart against the same $R_t = 1$ reference
line. Not executed in R; the printed range/robustness statement should
be checked once `estimate_r0()` can actually be run.

### Commit SHA / `renv.lock` / `sessionInfo()` pinning

Neither file previously contained any reproducibility-pinning code at
all (verified by search - no prior `renv`/SHA/`sessionInfo` references
existed anywhere in either file). Added to both:
- a check of the actually-installed `tsgc` version and (where
  available via `install_github()`/`install_git()`) commit/remote SHA,
  with an optional `EXPECTED_TSGC_COMMIT` to pin/verify against;
- a `sessionInfo()` capture written to disk;
- an `renv::snapshot()` call (guarded by `requireNamespace("renv")`)
  to pin dependency versions, not just `tsgc` itself.

In the script this is gated behind `SAVE_TABLES` (consistent with the
rest of the script's file-writing conventions); in the vignette it is
a separate chunk with `eval = FALSE` by default, since it writes to
disk and the reviewers' request was for the *capability* to pin a run,
not for every vignette render to write `sessionInfo.txt`/`renv.lock`
as a side effect. None of this was executed against the actually
installed package in this environment (no R available here) - only
the pinning/verification infrastructure was wired up; the printed
version/commit and the written `renv.lock`/`sessionInfo.txt` still
need to be generated/checked by actually running this code.

## Round 4 corrections (31 July 2026, same day)

Re-examined the three items previously left in "Outstanding" on the
grounds that they needed source files not in this batch or an R
environment. On rereading, `xts_to_idx()` and the leading-indicator
alignment code (`filterResultsLI.R`) were already in this batch - the
"not available" note in the round-1/round-3 changelog entries was
incorrect. All three are now genuinely resolved rather than deferred.

### `xts_to_idx()` duplicate-date behaviour - already resolved, test corrected

`xts_to_idx()`'s actual source (`utils.R`, lines 220-224) checks
`if (anyDuplicated(idx)) stop("xts_to_idx: x's index contains ",
"duplicate values.")` unconditionally, before any other logic runs.
A duplicated index date always errors clearly - this was never an
open defect. The round-1 test's hedged "error or de-duplicate" framing
and its claim that the source was unavailable were both wrong.
**Fix:** rewrote the test to assert the actual, confirmed behaviour
directly (`expect_error(xts_to_idx(malformed_xts), "duplicate
values")`), and corrected the surrounding comment.

### `SSModelLeadingIndicator`/`filterResultsLI.R` alignment logic - genuine defect found and fixed

Reading `filterResultsLI.R`'s `predict_all()` (not deferred this time)
found that `get_timeframe()` (`utils.R`) does not error on an
out-of-range request - it clamps: `start <- max(start, rng[1]); end <-
min(end, rng[2])`, silently returning fewer rows than requested rather
than failing. `predict_all()` calls
`get_timeframe(xpred_lead.new/xpred_targ.new, end+1, end+n.ahead)` and
then assigns the result positionally into a fixed-size
`n.ahead`-length array (`newZ[1, 1:d1, ] <- t(xpred_lead.new.subset)`)
with no check that the returned window actually has `n.ahead` rows. A
regressor series that runs out before the forecast horizon ends (a
genuinely realistic case - e.g. a weather feed that stops a few days
short) would be silently truncated and then recycled/misaligned into
that array, rather than raising the "fail clearly if any date is
missing" behaviour the critical assessments asked for.

**Fix (`filterResultsLI.R`):** both the `xpred_lead.new` and
`xpred_targ.new` branches of `predict_all()` now check that
`get_timeframe(..., end+1, end+n.ahead)` returned exactly `n.ahead`
positions before using the result, and `stop()` with a specific
message (including the expected/actual row counts and position range)
if not.

**Also discovered:** `idx_series` objects cannot represent an internal
gap at all - `idx_series(data, start)` always implies the contiguous
run of positions `start:(start+nrow(data)-1)`, and `idx_rbind()`
requires its two inputs to be adjoining and errors otherwise. The
round-1 "missing-date" test built its gappy series via
`idx_series(..., start = keep_pos[1])` after dropping a row, which
does not preserve the gap - it silently re-indexes every row after the
gap to close it up, shifting them to the wrong position/date instead
of leaving a hole. That test was not actually exercising a missing
date; it was, at best, exercising an off-by-one misalignment via an
indirect path. **Fix (`test-filterResultsLI.R`):** replaced it with a
test that truncates the regressor series a few positions before the
forecast horizon ends (the realistic failure mode, since a genuine
mid-series gap is impossible to construct), and added a second test
that pins `get_timeframe()`'s actual clamping behaviour directly, so
the reason `predict_all()` needs its own explicit length check (rather
than relying on `get_timeframe()` to error) is verified, not just
asserted in a comment.

None of this was executed in R (no environment available here); the
fix was derived by tracing the actual call graph
(`predict_all()` → `get_timeframe()` → its documented clamping
behaviour) and re-deriving the array-assignment consequence by hand,
not by running the code. A live test run is the way to confirm the
new `stop()` messages fire as written and that the existing
exact-date-alignment test (which supplies a full-length series and
still passes) is unaffected by the added length check.

## Outstanding (not addressed in this delivery)

Only the items genuinely outside this batch's scope remain:

- **Realised (oracle) future weather** — disclosed prominently (in
  comments and in every affected figure title) since round 1, but not
  resolved to an operational validation: that would require archived
  point/ensemble weather forecasts as of each origin date, which were
  never supplied in any batch and cannot be substituted with what's
  available here.
- **Full test suite / `R CMD check --as-cran` / clean-session render**
  — none of this was executed, since no R environment is available
  here. Every fix across all four rounds was verified by static
  inspection (reading source, tracing call graphs, and independent
  arithmetic re-derivation of window/threshold bounds), not by running
  the R code itself. This matters most for the round-4 fixes above,
  where the previous rounds' "not available"/"outstanding" framing was
  itself found to be incorrect on closer reading - an actual run
  (ideally a full `R CMD check --as-cran` plus a clean-session
  vignette render) is the only way to confirm the new `stop()`
  messages and rewritten tests behave as derived here.

## Files in this delivery (rounds 1-4, cumulative)

| File | Nature of changes |
|---|---|
| `tsgc_replication_script.R` | Core fixes (#1-2, 7-10), reinit fixes (#3-6), reporting/diagnostics infrastructure (#15-18), CV redesign + disjointness fix (#19-21), small-sample labelling (#22-23), gen_int sensitivity (#24), reproducibility pinning |
| `tsgc_vignette.Rmd` | Reinit fixes (#3-6), warning suppression (#15), CV redesign + disjointness fix (#19-21), gen_int sensitivity (#24), reproducibility pinning |
| `idx_calendar.R` | Date-class preservation fix (#11) |
| `test-idx_calendar.R` | Regression tests for #11 |
| `filterResults.R` | sMAPE docstring (#14) |
| `filterResultsLI.R` | sMAPE docstring (#14); **round 4:** explicit forecast-horizon-length check in `predict_all()` for both `xpred_lead.new` and `xpred_targ.new` |
| `test-filterResultsLI.R` | Regression tests for #13, driven through the model's real fields; numerical sMAPE unit test for #14; **round 4:** corrected duplicate-date test (pins `xts_to_idx()`'s actual unconditional error), replaced the flawed gappy-series test with a genuine short-series test, added a test pinning `get_timeframe()`'s clamping behaviour |

## Files in this delivery (round 1, for reference)

| File | Nature of changes |
|---|---|
| `tsgc_replication_script.R` | Core fixes (#1–2, 7–10), reinit fixes (#3–6), reporting/diagnostics infrastructure (#15–18), CV redesign (#19–21) |
| `tsgc_vignette.Rmd` | Reinit fixes (#3–6), warning suppression (#15), CV redesign (#19–21) |
| `idx_calendar.R` | Date-class preservation fix (#11) |
| `test-idx_calendar.R` | Regression tests for #11 |
| `filterResults.R` | sMAPE docstring (#14) |
| `filterResultsLI.R` | sMAPE docstring (#14) |
| `test-filterResultsLI.R` | Regression tests for #13, now driven through the model's real fields; numerical sMAPE unit test for #14 |