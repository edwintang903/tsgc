### ==========================================================================
### Smoke-test script for the idx_series refactor of tsgc - LEADING
### INDICATOR model (SSModelLeadingIndicator / FilterResultsLI)
###
### This is NOT a testthat suite. It is a plain script that exercises every
### piece of the new LI code, printing results at each step, so a human (or
### an LLM) can read the transcript and sanity-check behaviour by eye.
###
### Requires the same source directory as test_script.R (idx_series.R,
### utils.R, accessorFns.R) plus SSModelLeadingIndicator.R and
### filterResultsLI.R.
###
### Usage: Rscript test_script_LI.R 2>&1 | tee test_output_LI.txt
### ==========================================================================

section <- function(title) {
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat(title, "\n")
  cat(strrep("=", 80), "\n")
}

subsection <- function(title) {
  cat("\n")
  cat(strrep("-", 60), "\n")
  cat(title, "\n")
  cat(strrep("-", 60), "\n")
}

expect_true <- function(desc, cond) {
  status <- if (isTRUE(cond)) "PASS" else "**FAIL**"
  cat(sprintf("[%s] %s\n", status, desc))
}

expect_error <- function(desc, expr) {
  ok <- tryCatch({
    force(expr)
    FALSE
  }, error = function(e) {
    cat("    (caught expected error:", conditionMessage(e), ")\n")
    TRUE
  })
  status <- if (ok) "PASS" else "**FAIL**"
  cat(sprintf("[%s] %s (expected an error)\n", status, desc))
}

expect_equal <- function(desc, got, want, tol = 1e-8) {
  ok <- isTRUE(all.equal(got, want, tolerance = tol))
  status <- if (ok) "PASS" else "**FAIL**"
  cat(sprintf("[%s] %s\n", status, desc))
  if (!ok) {
    cat("    got: ", paste(capture.output(print(got)), collapse = "\n         "), "\n")
    cat("    want:", paste(capture.output(print(want)), collapse = "\n         "), "\n")
  }
}

### --------------------------------------------------------------------------
### 0. Load source
### --------------------------------------------------------------------------
section("0. LOADING SOURCE FILES")

pkgs <- c("KFAS", "magrittr", "methods", "abind", "purrr")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  package '%s' available: %s\n", p, ok))
}
# NB: SSModelLeadingIndicator.R/filterResultsLI.R call KFAS functions
# (SSModel, KFS, SSMtrend, SSMseasonal, SSMregression, fitSSM,
# predict.SSModel, ...) unqualified. Inside an installed package this
# works via NAMESPACE @importFrom, but since we're only source()-ing
# these files here, KFAS must be attached with library() for those bare
# calls to resolve. Same story for magrittr's %>% and purrr::partial.
library(KFAS)
library(magrittr)
library(methods)
library(abind)
library(purrr)

src_dir <- "./R"  # adjust if needed
src_files <- c("idx_series.R", "utils.R", "accessorFns.R",
               "SSModelLeadingIndicator.R", "filterResultsLI.R")
for (f in src_files) {
  cat("Sourcing", f, "...\n")
  source(file.path(src_dir, f))
}
cat("All source files loaded.\n")

subsection("0.1 Key objects exist after sourcing (catches silent load failures)")
expect_true("SSModelLeadingIndicator exists", exists("SSModelLeadingIndicator"))
expect_true("SSModelLeadingIndicator is a refObjectGenerator", inherits(SSModelLeadingIndicator, "refObjectGenerator"))
expect_true("FilterResultsLI exists", exists("FilterResultsLI"))
expect_true("FilterResultsLI is a refObjectGenerator", inherits(FilterResultsLI, "refObjectGenerator"))
if (!exists("SSModelLeadingIndicator") || !inherits(SSModelLeadingIndicator, "refObjectGenerator")) {
  cat("\n*** SSModelLeadingIndicator was not created correctly. Stopping here. ***\n")
  stop("SSModelLeadingIndicator not found after sourcing - see messages above.")
}

### --------------------------------------------------------------------------
### 1. Simulate data
### --------------------------------------------------------------------------
section("1. SIMULATING LEADING-INDICATOR DATA")

set.seed(42)
n_obs <- 150
true_lag <- 5

# Simulate a leading indicator series, then construct a target series whose
# growth rate loosely tracks the leading indicator's growth rate from
# `true_lag` positions earlier, so the model has genuine signal to recover.
lead_rate <- pmax(0.03 + 0.02 * sin(seq_len(n_obs) / 15) + rnorm(n_obs, 0, 0.01), 0.005)
lead_increments <- pmax(round(40 * cumprod(1 + lead_rate) * (1 + rnorm(n_obs, 0, 0.05))), 1)
cLead_full <- cumsum(lead_increments)

targ_rate <- c(rep(0.02, true_lag), lead_rate[1:(n_obs - true_lag)]) * 0.8 + 0.01
targ_increments <- pmax(round(25 * cumprod(1 + targ_rate) * (1 + rnorm(n_obs, 0, 0.05))), 1)
cTarg_full <- cumsum(targ_increments)

Y_full <- idx_series(cbind(cLead_full, cTarg_full), start = 1L)
cat("Simulated leading-indicator series Y_full, n =", length(Y_full), "\n")
cat("First 6 rows:\n")
print(head(idx_values(Y_full)))
cat("Last 6 rows:\n")
print(tail(idx_values(Y_full)))

### --------------------------------------------------------------------------
### 2. Basic model: fixed q, no seasonality
### --------------------------------------------------------------------------
section("2. Basic model: fixed q, no seasonality")

model_basic <- SSModelLeadingIndicator$new(
  Y = Y_full, n.lag = true_lag, sea.period = 0, q = 0.01,
  LeadIndCol = 1, start = 1, end = 120
)
cat("\n-- print(model_basic) --\n")
model_basic$print()
cat("\n-- summary(model_basic) --\n")
model_basic$summary()

res_basic <- model_basic$estimate()
cat("\nClass of res_basic:", class(res_basic), "\n")
expect_true("estimate() returns a FilterResultsLI", inherits(res_basic, "FilterResultsLI"))

cat("\n-- print(res_basic) --\n")
res_basic$print()
cat("\n-- summary(res_basic) --\n")
res_basic$summary()

cat("\nres_basic$start / end:", res_basic$start, "/", res_basic$end, "\n")
expect_true("res_basic$start is a sensible position", res_basic$start >= 1)
expect_true("res_basic$end matches model end (120)", res_basic$end == 120)

### --------------------------------------------------------------------------
### 3. predict_level / predict_all
### --------------------------------------------------------------------------
section("3. predict_level / predict_all")

subsection("3.1 predict_level (sea.on = FALSE)")
pl <- res_basic$predict_level(n.ahead = 10, confidence.level = 0.68, sea.on = FALSE)
cat("predict_level output (class:", class(pl), "):\n")
print(pl)
expect_true("predict_level returns idx_series", is_idx_series(pl))
expect_equal("predict_level ncol == 3 (forc/lwr/upr)", idx_ncol(pl), 3)
expect_equal("predict_level length == n.ahead", length(pl), 10)
pl_mat <- idx_values(pl)
expect_true("fit sits between lwr and upr for every forecast row",
            all(pl_mat[, "lwr"] <= pl_mat[, "forc"] & pl_mat[, "forc"] <= pl_mat[, "upr"]))
expect_true("no NA in predict_level output", !any(is.na(pl_mat)))

subsection("3.2 predict_level (n.ahead = 1, exercises the 'unity' branch)")
pl1 <- res_basic$predict_level(n.ahead = 1, sea.on = FALSE)
cat("predict_level(n.ahead=1) output:\n")
print(pl1)
expect_equal("predict_level(n.ahead=1) length == 1", length(pl1), 1)

subsection("3.3 predict_all (return.all = FALSE)")
pa <- res_basic$predict_all(n.ahead = 10, sea.on = FALSE, return.all = FALSE)
cat("predict_all names:", paste(names(pa), collapse = ", "), "\n")
cat("predict_all$y.hat:\n")
print(pa$y.hat)
cat("predict_all$level.t.t:\n")
print(pa$level.t.t)
cat("predict_all$slope.t.t:\n")
print(pa$slope.t.t)
expect_true("predict_all$y.hat is idx_series", is_idx_series(pa$y.hat))
expect_equal("predict_all$y.hat length == n.ahead", length(pa$y.hat), 10)

subsection("3.4 predict_all (return.all = TRUE)")
pa_full <- res_basic$predict_all(n.ahead = 10, sea.on = FALSE, return.all = TRUE)
expect_true("predict_all(return.all=TRUE) covers full sample + forecast",
            length(pa_full$y.hat) > length(pa$y.hat))

### --------------------------------------------------------------------------
### 4. get_growth_y / get_gy_ci
### --------------------------------------------------------------------------
section("4. get_growth_y / get_gy_ci")

gy <- res_basic$get_growth_y(smoothed = FALSE, return.components = FALSE)
cat("get_growth_y (filtered) head:\n")
print(head(idx_values(gy)))
expect_true("get_growth_y returns idx_series", is_idx_series(gy))

gy_comp <- res_basic$get_growth_y(smoothed = TRUE, return.components = TRUE)
cat("\nget_growth_y (smoothed, components) - list of length", length(gy_comp), "\n")
cat("gy.t head:\n"); print(head(idx_values(gy_comp[[1]])))
cat("g.t head:\n"); print(head(idx_values(gy_comp[[2]])))
cat("gamma.t head:\n"); print(head(idx_values(gy_comp[[3]])))

gy_ci <- res_basic$get_gy_ci(smoothed = FALSE, confidence.level = 0.68)
cat("\nget_gy_ci head:\n")
print(head(idx_values(gy_ci)))
expect_equal("get_gy_ci ncol == 3", idx_ncol(gy_ci), 3)
gy_ci_mat <- idx_values(gy_ci)
expect_true("get_gy_ci: fit sits between lower and upper (excluding degenerate diffuse-prior rows)",
            all(gy_ci_mat[-1, "lower"] <= gy_ci_mat[-1, "fit"] & gy_ci_mat[-1, "fit"] <= gy_ci_mat[-1, "upper"]))

### --------------------------------------------------------------------------
### 5. mapes
### --------------------------------------------------------------------------
section("5. mapes")

mp <- res_basic$mapes(n.ahead = 10, Y = Y_full)
cat("mapes() result:\n")
print(mp)
expect_true("mapes returns a list with mape/smape/mae/rmse/coverage",
            all(c("mape", "smape", "mae", "rmse", "coverage") %in% names(mp)))
expect_true("mapes$mape is a finite non-negative number", is.finite(mp$mape) && mp$mape >= 0)

### --------------------------------------------------------------------------
### 6. print_estimation_results (LaTeX table)
### --------------------------------------------------------------------------
section("6. print_estimation_results (LaTeX table)")

tryCatch({
  print(res_basic$print_estimation_results())
}, error = function(e) {
  cat("print_estimation_results() raised an error (likely missing 'kableExtra'):\n")
  cat("  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 7. Model with seasonality
### --------------------------------------------------------------------------
section("7. Model with seasonality")

model_sea <- SSModelLeadingIndicator$new(
  Y = Y_full, n.lag = true_lag, sea.period = 7, q = 0.01,
  LeadIndCol = 1, start = 1, end = 120
)
res_sea <- model_sea$estimate()
cat("res_sea class:", class(res_sea), "\n")
res_sea$print()

pl_sea <- res_sea$predict_level(n.ahead = 7, sea.on = TRUE)
cat("\npredict_level (seasonal model):\n")
print(pl_sea)
expect_true("seasonal model predict_level has no NA", !any(is.na(idx_values(pl_sea))))

tryCatch({
  print(res_sea$print_estimation_results())
}, error = function(e) {
  cat("print_estimation_results() [seasonal] raised an error:\n  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 8. Model with estimated q (q = NULL)
### --------------------------------------------------------------------------
section("8. Model with estimated q (q = NULL)")

model_estq <- SSModelLeadingIndicator$new(
  Y = Y_full, n.lag = true_lag, sea.period = 0, q = NULL,
  LeadIndCol = 1, start = 1, end = 100
)
res_estq <- model_estq$estimate()
cat("summary(res_estq):\n")
res_estq$summary()

### --------------------------------------------------------------------------
### 9. LeadIndCol = 2 (swap which column is the leading indicator)
### --------------------------------------------------------------------------
section("9. LeadIndCol = 2")

tryCatch({
  # Same data, but presented with target in column 1 and leading indicator
  # in column 2, using LeadIndCol=2 to tell the model which is which.
  Y_swapped <- idx_series(cbind(cTarg_full, cLead_full), start = 1L)
  model_swapped <- SSModelLeadingIndicator$new(
    Y = Y_swapped, n.lag = true_lag, sea.period = 0, q = 0.01,
    LeadIndCol = 2, start = 1, end = 120
  )
  res_swapped <- model_swapped$estimate()
  cat("res_swapped class:", class(res_swapped), "\n")
  res_swapped$print()
  pl_swapped <- res_swapped$predict_level(n.ahead = 5, sea.on = FALSE)
  cat("\npredict_level (LeadIndCol=2 model):\n")
  print(pl_swapped)
}, error = function(e) {
  cat("LeadIndCol=2 model raised an error:\n  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 10. Model with xpred_targ (exogenous predictor on target only)
### --------------------------------------------------------------------------
section("10. Model with xpred_targ (exogenous predictor on target)")

tryCatch({
  set.seed(11)
  xpred_targ_full <- idx_series(matrix(rnorm(n_obs), ncol = 1), start = 1L)
  model_xpred <- SSModelLeadingIndicator$new(
    Y = Y_full, n.lag = true_lag, sea.period = 0, q = 0.01,
    LeadIndCol = 1, xpred_targ = xpred_targ_full, start = 1, end = 100
  )
  res_xpred <- model_xpred$estimate()
  cat("res_xpred class:", class(res_xpred), "\n")
  res_xpred$print()

  res_xpred$xpred_targ.new <- get_timeframe(xpred_targ_full, 101, 150)
  pl_xpred <- res_xpred$predict_level(n.ahead = 10, sea.on = FALSE)
  cat("\npredict_level with xpred_targ:\n")
  print(pl_xpred)
}, error = function(e) {
  cat("xpred_targ model raised an error:\n  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 11. Model with xpred_lead (exogenous predictor on leading indicator)
### --------------------------------------------------------------------------
section("11. Model with xpred_lead (exogenous predictor on leading indicator)")

tryCatch({
  set.seed(13)
  xpred_lead_full <- idx_series(matrix(rnorm(n_obs), ncol = 1), start = 1L)
  model_xpred_lead <- SSModelLeadingIndicator$new(
    Y = Y_full, n.lag = true_lag, sea.period = 0, q = 0.01,
    LeadIndCol = 1, xpred_lead = xpred_lead_full, start = 1, end = 100
  )
  res_xpred_lead <- model_xpred_lead$estimate()
  cat("res_xpred_lead class:", class(res_xpred_lead), "\n")
  res_xpred_lead$print()

  res_xpred_lead$xpred_lead.new <- get_timeframe(xpred_lead_full, 101, 150)
  pl_xpred_lead <- res_xpred_lead$predict_level(n.ahead = 10, sea.on = FALSE)
  cat("\npredict_level with xpred_lead:\n")
  print(pl_xpred_lead)
}, error = function(e) {
  cat("xpred_lead model raised an error:\n  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 12. Error handling
### --------------------------------------------------------------------------
section("12. Error handling")

subsection("12.1 Invalid sea.period")
expect_error("sea.period == 1 should error",
             SSModelLeadingIndicator$new(Y = Y_full, n.lag = true_lag, sea.period = 1))
expect_error("sea.period < 0 should error",
             SSModelLeadingIndicator$new(Y = Y_full, n.lag = true_lag, sea.period = -1))

subsection("12.2 Invalid LeadIndCol")
expect_error("LeadIndCol not in {1,2} should error",
             SSModelLeadingIndicator$new(Y = Y_full, n.lag = true_lag, LeadIndCol = 3))

subsection("12.3 Invalid xpred type")
expect_error("xpred_lead must be idx_series or NULL",
             SSModelLeadingIndicator$new(Y = Y_full, n.lag = true_lag, xpred_lead = 1:10))
expect_error("xpred_targ must be idx_series or NULL",
             SSModelLeadingIndicator$new(Y = Y_full, n.lag = true_lag, xpred_targ = 1:10))

subsection("12.4 Non-increasing series after lagging")
bad_lead <- c(10, 20, 15, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150)
bad_targ <- cumsum(rep(5, length(bad_lead))) + 1
bad_Y <- idx_series(cbind(bad_lead, bad_targ), start = 1L)
expect_error("estimate() should error on non-increasing leading indicator",
             SSModelLeadingIndicator$new(Y = bad_Y, n.lag = 2, sea.period = 0)$estimate())

### --------------------------------------------------------------------------
### 13. cross_val (mixed model_list containing an LI model)
### --------------------------------------------------------------------------
section("13. cross_val with a leading-indicator model")

tryCatch({
  cv_models <- list(
    "LI_fixed_q" = SSModelLeadingIndicator$new(
      Y = Y_full, n.lag = true_lag, sea.period = 0, q = 0.01,
      LeadIndCol = 1, start = 1, end = 80),
    "LI_est_q" = SSModelLeadingIndicator$new(
      Y = Y_full, n.lag = true_lag, sea.period = 0, q = NULL,
      LeadIndCol = 1, start = 1, end = 80)
  )
  cv_results <- cross_val(
    Y = Y_full, model_list = cv_models, est.end = 80,
    n.ahead = 5, n.estimate = 2, gap = 5, criterion = "mape",
    LeadIndCol = 1
  )
  cat("cross_val results:\n")
  print(cv_results)
}, error = function(e) {
  cat("cross_val raised an error:\n  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 14. Accessor-style consistency checks
### --------------------------------------------------------------------------
section("14. ACCESSOR-STYLE CHECKS")

cat("output(res_basic) class:", class(output(res_basic)), "\n")
cat("modelKFS(output(res_basic)) class:", class(modelKFS(output(res_basic))), "\n")
cat("att(output(res_basic)) dim:", paste(dim(att(output(res_basic))), collapse = " x "), "\n")
cat("Ptt(output(res_basic)) dim:", paste(dim(Ptt(output(res_basic))), collapse = " x "), "\n")
cat("alphahat(output(res_basic)) dim:", paste(dim(alphahat(output(res_basic))), collapse = " x "), "\n")
cat("gety(modelKFS(output(res_basic))) dim:", paste(dim(gety(modelKFS(output(res_basic)))), collapse = " x "), "\n")
cat("seasonalComp(output(res_sea)):", !is.null(seasonalComp(output(res_sea))), "\n")

### --------------------------------------------------------------------------
### Done
### --------------------------------------------------------------------------
section("DONE")
cat("LI test script completed. Scroll up for any **FAIL** markers or uncaught errors.\n")
