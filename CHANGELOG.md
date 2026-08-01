# tsgc Replication — Correction Changelog

1 August 2026

## Corrections

1. **Gauteng `q` mismatch.** The forecast table was produced by a model
   with a freely estimated `q`, while the reported holdout accuracy came
   from a model with a fixed `q` — an apples-to-oranges comparison. Both
   now use the same fixed `q`, and all associated figures and tables were
   regenerated.
2. **Gauteng weather data overwritten by a rounded CSV.** A stray
   assignment copied 2-decimal-place CSV weather values over the
   full-precision series used for estimation, so every downstream figure
   silently used the rounded data. The assignment was removed and the
   CSV-derived series is now kept in its own object.
3. **Reinitialisation combined mismatched estimate/variance pairs.** The
   smoothed state estimate was paired with the filtered-state covariance
   rather than the matching smoothed covariance. Corrected to use one
   consistent pair throughout.
4. **Reinitialisation trigger compared the wrong quantities.** The
   crossing rule compared a lagged slope value against the threshold at
   the *current* date rather than each date's own threshold. Corrected
   so the signal and its threshold are always compared at the same date.
5. **Reinitialisation plot mislabelled its reference lines.** Lines shown
   as "standard-error bands around the slope" are k-SE *thresholds*, not
   confidence bands. Labels corrected.
6. **Reinitialisation overstated as a validated real-time rule.** Now
   explicitly disclosed as a retrospective, single-episode illustration.
   A smoothed estimate computed with data through the forecast origin
   cannot establish that the rule would have triggered in real time on
   an earlier date, since it incorporates information not yet available
   on that earlier date. The cross-validation redesign below addresses
   the closely related evaluation-design concern.
7. **Nintendo annual dates relabelled onto a fabricated calendar.** The
   selected quarterly observations were Q3 endpoints, not January
   observations, but had been assigned a made-up 1 January date. The
   series is now described as annual-frequency, Q3-endpoint cumulative
   sales, with its actual dates preserved.
8. **Filtered-estimate exports extended into the forecast period.**
   Files intended to hold only in-sample filtered estimates were not cut
   off at the end of the estimation sample. Now restricted accordingly.
9. **Diffuse-state rows exported with a false-certainty zero standard
   error.** Rows flagged as diffuse still reported a literal zero
   standard error underneath the flag, so any consumer not checking the
   flag would read zero uncertainty. Standard errors and interval bounds
   on flagged rows are now also set to missing.
10. **Multi-step calendars lost their Date/POSIXct class.** Positions
    converted back to dates under certain multi-step calendar patterns
    came back as raw numeric offsets rather than proper dates. Root
    cause and full resolution are described under Further fixes below,
    since the first attempted fix here was itself incomplete.
11. **No regression test for the multi-step calendar class-drop.**
    Added, alongside correcting three pre-existing tests that had
    asserted the wrong result type.
12. **England regressor-alignment tests didn't test the alignment
    logic.** The tests referenced a field that doesn't exist on the
    model object, so none of them actually exercised the real alignment
    path. Rewritten against the real fields, including a genuinely
    malformed input for the duplicate-date case.
13. **sMAPE formula and scale undocumented.** Documentation now states
    the exact formula and clarifies the implemented range is [0, 100],
    not the conventional [0, 200]. A numerical test pins the formula
    against a real fitted model's output, not just the formula in
    isolation.
14. **Global warning suppression hid real issues.** Removed, so warnings
    surface instead of being silently discarded.
15. **Two rendered figures were zero bytes despite a successful build.**
    Figure validation existed but was only run under a flag that
    defaults off. The save step now reports whether it actually wrote a
    file, and validation is tied to that outcome rather than to the
    flag.
16. **`summary()` reported only Length/Class/Mode.** Replaced with real
    diagnostics: estimated parameters, positions, and model states.
17. **Realised (oracle) future weather used without disclosure.** Now
    disclosed in comments and figure titles. Not resolved to an
    operational validation, which would require archived forecasts as of
    each origin date that were not supplied.
18. **An illustrative lag used downstream without explanation versus the
    cross-validation-selected lag.** Now disclosed as illustrative.
19. **Cross-validation selection and reporting blocks were not actually
    disjoint.** Forecast origins were assumed to step backward from the
    estimation end, but the underlying function steps forward, so the
    reporting fold's horizon could land on or past the selection block's
    boundary — a structural off-by-one, not a one-off coincidence.
    Corrected, with an explicit runtime check that the two blocks remain
    disjoint.
20. **No naive/random-walk benchmark in cross-validation.** Added.
21. **Annual leading-indicator accuracy based on only two observations.**
    Labelled explicitly as illustrative given the small sample.
22. **Interval coverage reported from very small evaluation samples.**
    Same treatment extended to all holdouts whose small size comes from
    each point being a whole year/quarter/month; daily-frequency
    holdouts were left as-is since that concern doesn't apply to them.
23. **No sensitivity analysis for the generation-interval assumption.**
    Added a sweep over a plausible range of values, reporting whether
    the qualitative conclusion is robust across it.
24. **No reproducibility pinning.** Added a check of the installed
    package version/commit against an expected value, a captured session
    summary, and a guarded environment-snapshot step.

## Further fixes

These weren't on the original list above; they turned up while making
the changes described there.

_Items below the "1 August 2026 (later pass)" marker were found in a
subsequent, independent review of the same replication script, vignette,
and package source (not the original correction pass above), and are
dated separately for that reason._

- **The first fix to the multi-step calendar date-class bug was itself
  wrong.** It passed a static read-through but failed the actual test
  suite. The real root cause was one level deeper: date-to-POSIXct
  promotion was being triggered by a property of the individual step
  rather than by the calendar's own settings, so plain date-only steps
  never got promoted and ended up mislabelled with both the wrong class
  and the wrong numeric scale. Fixed at the correct layer and confirmed
  by a full test run, including the tests that had previously failed.
- **A diagnostics helper errored on models with more than one
  observation-equation variance.** A boundary check assumed a single
  scalar variance value; leading-indicator models have two, and the
  check crashed on the resulting matrix. Fixed to handle any number of
  variance parameters, and the same underlying issue was found and fixed
  in an unrelated summary method that had been silently misprinting the
  same kind of matrix rather than erroring on it.
- **A single-row cross-validation window failed with a misleading
  low-level error instead of a clear one.** At the tightest lag setting,
  the first selection fold could leave exactly one usable row, which a
  matrix operation silently collapsed to a plain vector, surfacing a
  confusing error several layers downstream. Fixed at the source, with
  an explicit, clearly worded check in its place.
- **A duplicate-date test asserted hedged, incorrect behaviour.** The
  underlying function actually errors unconditionally on a duplicated
  date; this had never been an open defect. The test was rewritten to
  assert that directly.
- **A forecasting function didn't check regressor length before use.**
  An internal helper clamps rather than errors when asked for more data
  than is available, silently returning fewer rows. The function
  consuming its output assigned that result into a fixed-size array
  with no length check, so a regressor series running out before the
  forecast horizon ended would be silently truncated and misaligned.
  This was fixed and tested for the leading-indicator model
  (`xpred_lead.new`/`xpred_targ.new`), but the equivalent path for the
  plain Gompertz model (`xpred.new`) was missed in that pass: it still
  had no such check. Added the same explicit row-count check and error
  message there, plus a matching test, so both places now check the row
  count and fail clearly if it's short, as originally intended.
- **A "missing regressor date" test didn't actually test a missing
  date.** The data structure in question can't represent an internal
  gap at all, so the original approach of dropping a row silently
  re-indexed everything after it instead of leaving a hole. Replaced
  with a test of the realistic failure mode (a series that runs out
  before the forecast horizon ends), plus a second test pinning the
  clamping behaviour directly.
- **A diagnostics helper useful beyond this replication was promoted
  into the package itself**, rather than staying local to this script,
  along with the shared boundary-check logic it depends on. This also
  made both properly testable, closing the gap that a package-level
  test suite couldn't previously reach.
- **The vignette's weather-regressor example for England carried no
  oracle-forecast disclosure.** The same realised (not archived)
  future weather already disclosed for Gauteng and for England
  elsewhere in the replication script was used here too, but the
  vignette's prose, code comments, and figure title/caption were
  silent about it. Added the same caveat paragraph, code comment, and
  "(weather, oracle/realised)" figure labelling used elsewhere.
- **The vignette never demonstrated the corrected `summary()` output.**
  The reviewed defect was that all `summary()` calls in the vignette
  printed only Length/Class/Mode; the fix removed those calls rather
  than showing the corrected output, so a reader of the vignette alone
  saw no evidence `summary()` had been fixed. Added a "Model
  diagnostics" demonstration after both the baseline Gompertz model and
  the England leading-indicator model, calling `summary()` and
  `print_model_diagnostics()` on each - the latter specifically noting
  the matrix-valued `H` case for the leading-indicator model.
- **A vignette appendix silently plotted the wrong model.** The
  axis-mode appendix's helper function referenced the ambient `res`
  object and its caption claimed to reuse "the Gauteng fourteen-day
  new-cases forecast from Figure 3," but `res` had since been
  reassigned by the Reinitialization section (to the longer,
  25-June-end re-estimate) earlier in the same document. The appendix
  therefore rendered figures from the wrong model while describing them
  as the Figure 3 model. Fixed by re-fitting the original Figure 3
  specification under its own name (`res.axis.demo`) local to that
  appendix, rather than depending on what the shared `res` object
  happened to hold by that point in the document.
- **The reproducibility-pinning block referenced objects before they
  were defined.** The tsgc-commit/version check block called
  `ensure_dir(results_dir)` and `renv::snapshot(project = base_path)`,
  but `ensure_dir`, `results_dir`, and `base_path` were all defined
  later in the same script (Section 1.4). With `SAVE_TABLES = TRUE`
  this would fail immediately with an object/function-not-found error,
  so the reproducibility feature credited elsewhere in this changelog
  could not actually run as originally ordered. Moved the block to
  after Section 1.4, where all three are defined.
- **Leading-indicator forecast intervals now honour the requested 
  confidence level in all cases.** A code path in 
  `FilterResultsLI$predict_level()` ignored the user-supplied 
  `confidence.level` when seasonal adjustment was enabled, silently 
  falling back to the default interval width. Fixed so the requested 
  confidence level is applied consistently.
- **`plot_holdout()` now applies the supplied plot caption.** The 
  `caption` argument was accepted but never passed through to the 
  underlying `ggplot2` call, so user-specified captions were silently 
  ignored.
- **`SSModelLeadingIndicator()` now validates `n.lag` on construction.** 
  Invalid lag specifications (for example negative or non-integer values) 
  now fail immediately with a clear error rather than surfacing later 
  during estimation.

- **Corrected a non-runnable vignette example.** The weather-regressor 
  CSV example referenced an undefined `res_weather` object; the missing 
  model-fitting step has been restored so the example can be run as 
  presented.

## Out of scope

- **Annual leading-indicator accuracy from two observations** and
  **interval coverage from very small evaluation samples** — labelling
  fixes are above; the underlying small-sample limitation itself isn't
  fixable without more data.