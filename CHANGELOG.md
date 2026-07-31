# tsgc Replication — Correction Changelog

31 July 2026

Scope: `tsgc_replication_script.R`, `tsgc_vignette.Rmd`, `idx_calendar.R`,
`test-idx_calendar.R`, `filterResultsLI.R`, `test-filterResultsLI.R`,
`filterResults.R`, `SSModelLeadingIndicator.R`. Verified by direct
inspection of source; where noted, also verified by an actual R run.

## Fixed

1. **Gauteng `q` mismatch** — the forecast used a freely estimated `q`
   while the holdout evaluation used a fixed `q`, giving an
   apples-to-oranges accuracy comparison. *(`tsgc_replication_script.R`)*
2. **Gauteng weather CSV overwrote the full-precision data** — even
   after an initial fix kept the CSV-derived series in its own object,
   a leftover assignment line still copied the rounded (2 d.p.)
   CSV values over the full-precision series, so every downstream
   figure used the rounded data regardless. The stray assignment was
   removed. *(`tsgc_replication_script.R`)*
3. **Reinitialisation paired smoothed estimates with filtered variance**
   — `alphahat` (smoothed) was combined with the filtered `P`/`Ptt`
   instead of the matching smoothed covariance. *(`tsgc_replication_script.R`,
   `tsgc_vignette.Rmd`)*
4. **Reinitialisation trigger compared the wrong quantities** — the
   crossing rule compared a lagged slope value to the *current*
   threshold rather than comparing the signal to its own threshold at
   the same date. *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*
5. **Reinitialisation plot mislabelled its reference lines** — lines
   were labelled "SE bands around the slope" when they are k-SE
   *thresholds*, not confidence bands. *(`tsgc_replication_script.R`,
   `tsgc_vignette.Rmd`)*
6. **Reinitialisation exercise overstated as a real-time trigger** — now
   disclosed as a retrospective, single-episode illustration rather
   than a validated real-time rule; the CV redesign below (#19–21)
   is a related fix. *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*
7. **Nintendo annual dates relabelled onto a fabricated calendar** — the
   series had been given a made-up 1 Jan date rather than its actual
   annual dates. *(`tsgc_replication_script.R`)*
8. **`*_filtered.csv` files silently extended into the forecast period**
   — filtered-estimate exports weren't cut off at the end of the
   estimation sample. *(`tsgc_replication_script.R`)*
9. **Diffuse-state rows exported with a false-certainty zero SE** — an
    initial fix added a `diffuse_flag` column but still exported the
    literal zero standard error underneath it, so any consumer that
    didn't check the flag would still read zero uncertainty.
    `write_idx_csv()` now also sets the SE and any `_lower_`/`_upper_`
    bound columns to `NA` on flagged rows. *(`tsgc_replication_script.R`)*
10. **`idx_to_date()` dropped Date/POSIXct class for multi-step
    calendars** — the real bug was one level deeper than an initial
    attribute-copying fix assumed (see note below).
    *(`idx_calendar.R`)*
11. **No regression test for the multi-step calendar class-drop bug** —
    added; also corrected three pre-existing, unrelated tests that
    asserted a plain `Date` result on a `posixct = TRUE` calendar.
    *(`test-idx_calendar.R`)*
12. **England regressor alignment tests didn't test the alignment logic**
    — the three tests assigned expected values to `res$xpred.new`, a
    field that doesn't exist on the model object (the real fields are
    `xpred_lead.new`/`xpred_targ.new`), so none of them actually
    exercised `predict_all()`'s alignment path. All three rewritten to
    use the real fields and, for the duplicate-date case, a genuinely
    malformed `xts` input. *(`test-filterResultsLI.R`)*
13. **sMAPE formula/scale undocumented** — docstrings now state the
    exact formula and clarify the implemented range is [0, 100], not
    the conventional [0, 200]; added a numerical unit test pinning the
    formula directly. *(`filterResults.R`, `filterResultsLI.R`,
    `test-filterResultsLI.R`)*
14. **Global warning suppression hid real issues** — `warning = FALSE`
    removed so warnings surface instead of being silently discarded.
    *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*
15. **Two rendered figures were zero bytes despite a successful build**
    — `validate_saved_figure()` was correctly implemented but only ever
    called inside `if (SAVE_PLOTS) ...`, and `SAVE_PLOTS` defaults to
    `FALSE`, so the check never ran by default. `safe_ggsave()` now
    returns whether it actually wrote a file, and validation is tied to
    that return value; a comment now flags that `SAVE_PLOTS <- TRUE` is
    required to actually exercise render validation.
    *(`tsgc_replication_script.R`)*
16. **`summary()` reported only Length/Class/Mode** — replaced with real
    diagnostics. *(`tsgc_replication_script.R`)*
17. **Realised (oracle) future weather used without disclosure** — now
    disclosed in comments and figure titles. Not resolved to an
    operational validation — that would need archived point/ensemble
    forecasts as of each origin date, which weren't supplied in any
    batch. *(`tsgc_replication_script.R`)*
18. **Lag 14 used downstream with no explanation vs. the CV-selected
    lag** — now disclosed as illustrative. *(`tsgc_replication_script.R`,
    `tsgc_vignette.Rmd`)*
19. **Cross-validation SELECTION and REPORTING blocks were not actually
    disjoint** — origins were assumed to step *backward* from
    `est.end`, but `cross_val()` actually steps forward. This made the
    REPORTING fold's forecast horizon land exactly on (or past) the
    SELECTION block's boundary rather than strictly before it — a
    structural off-by-one, not a one-off coincidence. `select.end.cv`
    now subtracts `n.select.cv * gap.cv` instead of a single `gap.cv`,
    and an explicit `stopifnot()` checks disjointness at runtime.
    *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*
20. **No naive/random-walk benchmark in cross-validation** — added.
    *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*
21. **Annual leading-indicator accuracy based on only two observations**
    — added an explicit small-sample caveat and relabelled the plot
    title "(illustrative only – n=2 holdout observations)".
    *(`tsgc_replication_script.R`)*
22. **Interval coverage reported from very small evaluation samples** —
    same treatment extended to the quarterly Wii (n=4), Wii→Switch
    (n=8), monthly Plus500 (n=4), and DEGIRO→AvaTrade (n=4) holdouts.
    Daily-frequency (7- and 14-day) holdouts were left as-is, since the
    concern was specifically about samples that are small because each
    point is a whole year/quarter/month. *(`tsgc_replication_script.R`)*
23. **No sensitivity analysis for the generation-interval assumption** —
    added a `gen_int` sweep over `{3,4,5,6,7}` days, reporting the range
    of $R_t$ estimates and whether the above/below-1 conclusion is
    robust across the grid, plus a combined chart.
    *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*
24. **No reproducibility pinning (commit SHA / `renv.lock` /
    `sessionInfo()`)** — added a check of the installed `tsgc`
    version/commit against an optional `EXPECTED_TSGC_COMMIT`, a
    `sessionInfo()` capture, and a guarded `renv::snapshot()` call.
    Gated behind `SAVE_TABLES` in the script; a separate `eval = FALSE`
    chunk in the vignette (writes to disk, so not run on every render
    by default). *(`tsgc_replication_script.R`, `tsgc_vignette.Rmd`)*

## Fixed — found during a later re-audit, not from the original review list

25. **`idx_to_date()` fix above was itself wrong** — the initial fix
    (`attributes(out) <- attributes(cal$anchor)`) passed static review
    but failed the actual test suite (6 failures: `POSIXct` check
    returning `FALSE`, plus a `strptime`/invalid-`tz` error). Root
    cause: `idx_step_add()` only promotes `Date` to `POSIXct` when the
    *step itself* has a nonzero sub-day component — it never checks the
    calendar's `posixct` flag. Every failing test used date-only steps,
    so promotion never happened, and reused `Date` attributes
    mislabelled the result with both the wrong class and the wrong
    numeric scale (seconds-since-epoch read as days-since-epoch).
    `idx_multi_step_offset_to_date()` now promotes the anchor to
    `POSIXct` up front based on `cal$posixct` itself. Confirmed by an
    actual test run — the 6 failing tests plus 3 related corrected
    tests all pass. *(`idx_calendar.R`)*
26. **`print_model_diagnostics()` errored on multi-parameter `H`** — a
    live render surfaced `'length = 4' in coercion to 'logical(1)'`
    because the leading-indicator model's `H` is a 2×2 matrix, not a
    scalar, and the boundary check used `&&`. Replaced
    `!is.na(H) && H <= boundary_tol` with
    `!anyNA(H) && any(H <= boundary_tol)`. Confirmed by an actual run
    that previously errored at this line and now completes.
27. **Single-row CV estimation window silently malformed, surfaced as an
    opaque KFAS error** — an actual `cross_val()` run over
    `Lag1..Lag21` failed deep inside `SSModel()` with a misleading
    "Misspecified H" error. Root cause: at `n.lag = 21`, the first
    SELECTION fold leaves exactly one usable row, and a matrix subset
    without `drop = FALSE` silently collapsed it to a plain vector.
    Added `drop = FALSE`, plus an explicit `if (nrow(data_mat) < 2)
    stop(...)` naming `n.lag`, the row count, and the remedy.
    *(`SSModelLeadingIndicator.R`)*
28. **`xts_to_idx()` duplicate-date test was hedged/wrong** — the
    function's real source unconditionally errors on a duplicated
    index; this was never actually an open defect. The test was
    rewritten to assert that behaviour directly instead of hedging
    between "errors or de-duplicates." *(`test-filterResultsLI.R`)*
29. **`predict_all()` didn't check regressor length before array
    assignment** — `get_timeframe()` clamps rather than errors on an
    out-of-range request, silently returning fewer rows than asked
    for. `predict_all()` assigned that result positionally into a
    fixed-size `n.ahead`-length array with no length check, so a
    regressor series that runs out before the forecast horizon ends
    would be silently truncated and misaligned. Both
    `xpred_lead.new`/`xpred_targ.new` branches now check the returned
    row count and `stop()` with the expected/actual counts if it's
    short. *(`filterResultsLI.R`)*
30. **Round-1 "missing regressor date" test didn't test a missing date**
    — `idx_series` objects can't represent an internal gap at all
    (`idx_series(data, start)` always implies a contiguous run of
    positions), so the original test's approach of dropping a row
    silently re-indexed everything after it instead of leaving a hole.
    Replaced with a test that truncates the series before the forecast
    horizon ends (the realistic failure mode), and added a second test
    pinning `get_timeframe()`'s actual clamping behaviour.
    *(`test-filterResultsLI.R`)*

## Out of scope

- **Annual leading-indicator accuracy from two observations** and
  **interval coverage from very small evaluation samples** —
  labelling fixes are above (#22–23); the underlying small-sample
  limitation itself isn't fixable without more data.