### ==========================================================================
### Smoke-test script for the idx_series refactor of tsgc
###
### This is NOT a testthat suite. It is a plain script that exercises every
### piece of the new code, printing results at each step, so a human (or an
### LLM) can read the transcript and sanity-check behaviour by eye.
###
### Usage: Rscript test_script.R 2>&1 | tee test_output.txt
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
# NB: SSModelDynGompertz.R/filterResults.R call KFAS functions (SSModel,
# KFS, SSMtrend, SSMregression, fitSSM, predict.SSModel, ...) unqualified.
# Inside an installed package this works via NAMESPACE @importFrom, but
# since we're only source()-ing these files here, KFAS must be attached
# with library() for those bare calls to resolve. Same story for
# magrittr's %>% and abind::abind (abind is already namespace-qualified
# in the source, so it doesn't strictly need attaching, but there's no
# harm in it being attached too).
library(KFAS)
library(magrittr)
library(methods)
library(abind)

src_dir <- "./R"  # adjust if needed
src_files <- c("idx_series.R", "utils.R", "SSModelDynGompertz.R", "filterResults.R", "accessorFns.R")
for (f in src_files) {
  cat("Sourcing", f, "...\n")
  source(file.path(src_dir, f))
}
cat("All source files loaded.\n")

subsection("0.1 S3 dispatch sanity check (post source())")
.tmp <- idx_series(c(1, 2, 3), start = 1L)
expect_true("as.numeric.idx_series is registered", !is.null(getS3method("as.numeric", "idx_series", optional = TRUE)))
expect_true("as.matrix.idx_series is registered", !is.null(getS3method("as.matrix", "idx_series", optional = TRUE)))
expect_true("print.idx_series is registered", !is.null(getS3method("print", "idx_series", optional = TRUE)))
rm(.tmp)

subsection("0.2 Key objects exist after sourcing (catches silent load failures)")
expect_true("SSModelDynamicGompertz exists", exists("SSModelDynamicGompertz"))
expect_true("SSModelDynamicGompertz is a refObjectGenerator", inherits(SSModelDynamicGompertz, "refObjectGenerator"))
expect_true("FilterResults exists", exists("FilterResults"))
expect_true("FilterResults is a refObjectGenerator", inherits(FilterResults, "refObjectGenerator"))
if (!exists("SSModelDynamicGompertz") || !inherits(SSModelDynamicGompertz, "refObjectGenerator")) {
  cat("\n*** SSModelDynamicGompertz was not created correctly. Common causes: ***\n")
  cat("  - setRefClass() call in SSModelDynGompertz.R errored/warned during sourcing\n")
  cat("    (check for messages printed just above the 'Sourcing SSModelDynGompertz.R ...' line)\n")
  cat("  - the file was sourced into a different environment than globalenv()\n")
  cat("  - a leftover browser()/debug session is intercepting execution\n")
  cat("Stopping here so the rest of the script doesn't cascade-fail.\n")
  stop("SSModelDynamicGompertz not found after sourcing - see diagnostics above.")
}

### --------------------------------------------------------------------------
### 1. idx_series basics
### --------------------------------------------------------------------------
section("1. idx_series BASICS")

subsection("1.1 Construction (vector)")
x <- idx_series(c(10, 20, 35, 60, 90), start = 5L)
print(x)
expect_equal("length(x) == 5", length(x), 5)
expect_equal("idx_range(x) == c(5,9)", idx_range(x), c(5L, 9L))
expect_equal("idx_positions(x) == 5:9", idx_positions(x), 5:9)
expect_equal("idx_ncol(x) == 1", idx_ncol(x), 1)
expect_equal("idx_values(x) == data", idx_values(x), c(10, 20, 35, 60, 90))
expect_equal("as.numeric(x) == data (S3 dispatch check)", as.numeric(x), c(10, 20, 35, 60, 90))
expect_equal("as.double(x) == data (S3 dispatch check)", as.double(x), c(10, 20, 35, 60, 90))
expect_true("as.numeric(x) does NOT just return the raw list",
            !is.list(as.numeric(x)))

subsection("1.2 Construction (matrix)")
m <- idx_series(matrix(1:10, ncol = 2), start = 100L)
print(m)
expect_equal("length(m) == 5 (nrow)", length(m), 5)
expect_equal("idx_ncol(m) == 2", idx_ncol(m), 2)
expect_equal("idx_range(m) == c(100,104)", idx_range(m), c(100L, 104L))

subsection("1.3 Subsetting by position")
sub <- x[6:8]
print(sub)
expect_equal("x[6:8] start == 6", sub$start, 6L)
expect_equal("x[6:8] values", as.numeric(sub), c(20, 35, 60))

expect_error("subsetting out-of-range positions should error", x[1:3])

subsection("1.4 idx_diff / idx_lag")
d <- idx_diff(x, 1L)
print(d)
expect_equal("idx_diff(x,1) start == 6", d$start, 6L)
expect_equal("idx_diff(x,1) values == diffs", as.numeric(d), c(10, 15, 25, 30))

lg <- idx_lag(x, 1L)
print(lg)
expect_equal("idx_lag(x,1) start == 6", lg$start, 6L)
expect_equal("idx_lag(x,1)[6] == x[5]", as.numeric(lg[6:9]), as.numeric(x[5:8]))

subsection("1.5 idx_rbind / idx_cbind")
tail_part <- idx_series(c(150, 200), start = 10L)
combined <- idx_rbind(x, tail_part)
print(combined)
expect_equal("idx_rbind length", length(combined), 7)
expect_equal("idx_rbind start", combined$start, 5L)

y1 <- idx_series(c(1, 2, 3), start = 1L)
y2 <- idx_series(c(4, 5, 6), start = 1L)
cb <- idx_cbind(y1, y2)
print(cb)
expect_equal("idx_cbind ncol == 2", idx_ncol(cb), 2)

expect_error("idx_rbind with wrong start should error",
             idx_rbind(x, idx_series(c(1, 2), start = 999L)))

subsection("1.6 is_idx_series / as_idx_series")
expect_true("is_idx_series(x) TRUE", is_idx_series(x))
expect_true("is_idx_series(1:5) FALSE", !is_idx_series(1:5))
coerced <- as_idx_series(1:5, start = 3L)
expect_true("as_idx_series coerces plain vector", is_idx_series(coerced))
expect_equal("as_idx_series preserves start", coerced$start, 3L)
passthrough <- as_idx_series(x)
expect_true("as_idx_series passes through idx_series unchanged", identical(passthrough, x))

### --------------------------------------------------------------------------
### 2. utils.R
### --------------------------------------------------------------------------
section("2. utils.R HELPERS")

subsection("2.1 get_timeframe")
long_series <- idx_series(cumsum(rep(5, 30)) + 1, start = 1L)
gt <- get_timeframe(long_series, 5, 10)
print(gt)
expect_equal("get_timeframe start", gt$start, 5L)
expect_equal("get_timeframe length", length(gt), 6)

gt2 <- get_timeframe(long_series, 25)
expect_equal("get_timeframe with no end goes to series end", idx_range(gt2)[2], idx_range(long_series)[2])

expect_true("get_timeframe(NULL, ...) returns NULL", is.null(get_timeframe(NULL, 1, 5)))

subsection("2.2 df2ldl")
set.seed(42)
cum_series <- idx_series(cumsum(rpois(20, 8)) + 1, start = 1L)
ldl <- df2ldl(cum_series)
print(ldl)
expect_equal("df2ldl start", ldl$start, 2L)
expect_equal("df2ldl length", length(ldl), 19)
# manual check on the first value
manual_first <- log((as.numeric(cum_series[2]) - as.numeric(cum_series[1])) / as.numeric(cum_series[1]))
expect_equal("df2ldl first value matches manual calc", as.numeric(ldl[2]), manual_first)

subsection("2.3 reinitialise_dataframe")
reinit <- reinitialise_dataframe(cum_series, 10)
print(reinit)
expect_equal("reinitialise_dataframe start", reinit$start, 10L)
expect_equal("reinitialise_dataframe first value is the increment", as.numeric(reinit[10]),
             as.numeric(cum_series[10]) - as.numeric(cum_series[9]))

subsection("2.4 add_daily_ldl")
lead_targ <- idx_series(cbind(cumsum(rpois(20, 5)) + 1, cumsum(rpois(20, 9)) + 1), start = 1L)
ldl_list <- add_daily_ldl(lead_targ, LeadIndCol = 1)
cat("Names:", paste(names(ldl_list), collapse = ", "), "\n")
print(ldl_list$LDLlead)
print(ldl_list$LDLtarg)
expect_equal("add_daily_ldl LDLlead start", ldl_list$LDLlead$start, 2L)

subsection("2.5 argmax")
am <- argmax(cum_series)
print(am)
expect_true("argmax returns idx_series of length 1", is_idx_series(am) && length(am) == 1)

subsection("2.6 error handling")
expect_error("df2ldl on non-idx_series should error", df2ldl(1:10))
expect_error("df2ldl on 2-column series should error", df2ldl(lead_targ))
expect_error("reinitialise_dataframe with out-of-range idx should error",
             reinitialise_dataframe(cum_series, 100))

### --------------------------------------------------------------------------
### 3. SSModelDynamicGompertz + FilterResults
### --------------------------------------------------------------------------
section("3. SSModelDynamicGompertz / FilterResults")

set.seed(123)
n_obs <- 150
# Simulate a smooth-ish cumulative growth curve with a bit of noise, so
# increments stay strictly positive (a requirement of the model).
true_rate <- 0.06 * exp(-seq_len(n_obs) / 60) + 0.005
increments <- pmax(round(50 * cumprod(1 + true_rate) * (1 + rnorm(n_obs, 0, 0.05))), 1)
Y_full <- idx_series(cumsum(increments), start = 1L)
cat("Simulated cumulative series Y_full, n =", length(Y_full), "\n")
cat("First 10 values: ", paste(round(as.numeric(Y_full[1:10])), collapse = ", "), "\n")
cat("Last 5 values:   ", paste(round(as.numeric(Y_full[146:150])), collapse = ", "), "\n")

subsection("3.1 Basic model: fixed q, no seasonality, no ar1")
model_basic <- SSModelDynamicGompertz$new(Y = Y_full, q = 0.01, sea.period = 0, end = 120)
cat("\n-- print(model_basic) --\n")
model_basic$print()
cat("\n-- summary(model_basic) --\n")
model_basic$summary()

res_basic <- model_basic$estimate()
cat("\nClass of res_basic:", class(res_basic), "\n")
expect_true("estimate() returns a FilterResults", inherits(res_basic, "FilterResults"))

cat("\n-- print(res_basic) --\n")
res_basic$print()
cat("\n-- summary(res_basic) --\n")
res_basic$summary()

subsection("3.2 predict_level / predict_all")
pl <- res_basic$predict_level(n.ahead = 10, confidence.level = 0.68, sea.on = FALSE)
cat("\npredict_level output (class:", class(pl), "):\n")
print(pl)
expect_true("predict_level returns idx_series", is_idx_series(pl))
expect_equal("predict_level ncol == 3 (fit/lower/upper)", idx_ncol(pl), 3)
expect_equal("predict_level length == n.ahead", length(pl), 10)

pa <- res_basic$predict_all(n.ahead = 10, sea.on = FALSE, return.all = FALSE)
cat("\npredict_all names:", paste(names(pa), collapse = ", "), "\n")
cat("predict_all$y.hat:\n")
print(pa$y.hat)
cat("predict_all$level.t.t:\n")
print(pa$level.t.t)
cat("predict_all$slope.t.t:\n")
print(pa$slope.t.t)

pa_full <- res_basic$predict_all(n.ahead = 10, sea.on = FALSE, return.all = TRUE)
expect_true("predict_all(return.all=TRUE) covers full sample + forecast",
            length(pa_full$y.hat) > length(pa$y.hat))

subsection("3.3 get_growth_y / get_gy_ci")
gy <- res_basic$get_growth_y(smoothed = FALSE, return.components = FALSE)
cat("\nget_growth_y (filtered):\n")
print(gy)
expect_true("get_growth_y returns idx_series", is_idx_series(gy))

gy_comp <- res_basic$get_growth_y(smoothed = TRUE, return.components = TRUE)
cat("\nget_growth_y (smoothed, components) - list of length", length(gy_comp), "\n")
cat("gy.t head:\n"); print(head(as.numeric(gy_comp[[1]])))
cat("g.t head:\n"); print(head(as.numeric(gy_comp[[2]])))
cat("gamma.t head:\n"); print(head(as.numeric(gy_comp[[3]])))

gy_ci <- res_basic$get_gy_ci(smoothed = FALSE, confidence.level = 0.68)
cat("\nget_gy_ci output:\n")
print(gy_ci)
expect_equal("get_gy_ci ncol == 3", idx_ncol(gy_ci), 3)

subsection("3.4 mapes")
mp <- res_basic$mapes(n.ahead = 10, Y = Y_full)
cat("\nmapes() result:\n")
print(mp)
expect_true("mapes returns a list with mape/smape/mae/rmse/coverage",
            all(c("mape", "smape", "mae", "rmse", "coverage") %in% names(mp)))

subsection("3.5 print_estimation_results (LaTeX table)")
tryCatch({
  print(res_basic$print_estimation_results())
}, error = function(e) {
  cat("print_estimation_results() raised an error (likely missing 'kableExtra'):\n")
  cat("  ", conditionMessage(e), "\n")
})

subsection("3.6 Model with seasonality")
model_sea <- SSModelDynamicGompertz$new(Y = Y_full, q = 0.01, sea.period = 7, end = 120)
res_sea <- model_sea$estimate()
cat("\nres_sea class:", class(res_sea), "\n")
res_sea$print()
pl_sea <- res_sea$predict_level(n.ahead = 7, sea.on = TRUE)
cat("\npredict_level (seasonal model):\n")
print(pl_sea)

subsection("3.7 Model with estimated q (q = NULL)")
model_estq <- SSModelDynamicGompertz$new(Y = Y_full, q = NULL, sea.period = 0, end = 100)
res_estq <- model_estq$estimate()
cat("\nsummary(res_estq):\n")
res_estq$summary()

subsection("3.8 Model with ar1 = TRUE")
tryCatch({
  model_ar1 <- SSModelDynamicGompertz$new(Y = Y_full, q = 0.01, sea.period = 0, ar1 = TRUE, end = 100)
  res_ar1 <- model_ar1$estimate()
  cat("\nres_ar1 class:", class(res_ar1), "\n")
  res_ar1$print()
}, error = function(e) {
  cat("ar1 model raised an error:\n  ", conditionMessage(e), "\n")
})

subsection("3.9 Model with xpred (exogenous predictors)")
tryCatch({
  set.seed(7)
  xpred_mat <- matrix(rnorm(n_obs), ncol = 1)
  colnames(xpred_mat) <- "x1"
  xpred_full <- idx_series(xpred_mat, start = 1L)
  model_xpred <- SSModelDynamicGompertz$new(
    Y = Y_full, q = 0.01, sea.period = 0, xpred = xpred_full, end = 100
  )
  res_xpred <- model_xpred$estimate()
  cat("\nres_xpred class:", class(res_xpred), "\n")
  res_xpred$print()
  
  res_xpred$xpred.new <- get_timeframe(xpred_full, 101, 150)
  pl_xpred <- res_xpred$predict_level(n.ahead = 10, sea.on = FALSE)
  cat("\npredict_level with xpred:\n")
  print(pl_xpred)
}, error = function(e) {
  cat("xpred model raised an error:\n  ", conditionMessage(e), "\n")
})

subsection("3.10 Model with reinit.idx (reinitialisation)")
tryCatch({
  model_reinit <- SSModelDynamicGompertz$new(
    Y = Y_full, q = 0.01, sea.period = 0, reinit.idx = 60, end = 120
  )
  res_reinit <- model_reinit$estimate()
  cat("\nres_reinit class:", class(res_reinit), "\n")
  cat("res_reinit$index range:", range(res_reinit$index), "\n")
  res_reinit$print()
  res_reinit$summary()
  
  pl_reinit <- res_reinit$predict_level(n.ahead = 10, sea.on = FALSE)
  cat("\npredict_level (reinitialised model):\n")
  print(pl_reinit)
}, error = function(e) {
  cat("reinit model raised an error:\n  ", conditionMessage(e), "\n")
})

subsection("3.11 Error handling: non-increasing Y")
bad_Y <- idx_series(c(10, 20, 15, 40, 50, 60, 70, 80, 90, 100), start = 1L)
expect_error("estimate() should error on non-increasing Y",
             SSModelDynamicGompertz$new(Y = bad_Y, q = 0.01, sea.period = 0)$estimate())

subsection("3.12 Error handling: invalid sea.period")
expect_error("sea.period == 1 should error",
             SSModelDynamicGompertz$new(Y = Y_full, sea.period = 1))
expect_error("sea.period < 0 should error",
             SSModelDynamicGompertz$new(Y = Y_full, sea.period = -1))

### --------------------------------------------------------------------------
### 4. cross_val
### --------------------------------------------------------------------------
section("4. cross_val")

tryCatch({
  cv_models <- list(
    "q_fixed" = SSModelDynamicGompertz$new(Y = Y_full, q = 0.01, sea.period = 0, start = 1, end = 80),
    "q_est"   = SSModelDynamicGompertz$new(Y = Y_full, q = NULL, sea.period = 0, start = 1, end = 80)
  )
  cv_results <- cross_val(
    Y = Y_full, model_list = cv_models, est.end = 80,
    n.ahead = 5, n.estimate = 3, gap = 5, criterion = "mape"
  )
  cat("\ncross_val results:\n")
  print(cv_results)
}, error = function(e) {
  cat("cross_val raised an error:\n  ", conditionMessage(e), "\n")
})

### --------------------------------------------------------------------------
### 5. Accessor-style consistency checks (mirroring accessorFns.R usage)
### --------------------------------------------------------------------------
section("5. ACCESSOR-STYLE CHECKS")

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
cat("Test script completed. Scroll up for any **FAIL** markers or uncaught errors.\n")