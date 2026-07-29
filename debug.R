# --- Standalone debug script for MAPE/sMAPE/MAE/RMSE = NaN in Figures 4/8/9 ---
# Run this from the package root (devtools::load_all) or after library(tsgc).

library(tsgc)
library(dplyr)

# 1. Load data exactly as the vignette does -----------------------------
data(gauteng, package = "tsgc")
conv <- xts_to_idx(gauteng)
gauteng_idx <- conv$series
gauteng_cal <- conv$calendar
gauteng_cal$anchor_name <- "first recorded case"

Y   <- gauteng_idx
cal <- gauteng_cal

estimation.pos.start <- idx_to_pos(cal, "2021-02-01")
estimation.pos.end   <- idx_to_pos(cal, "2021-05-03")
n.forecasts <- 14
q <- 0.005
confidence.level <- 0.68

# 2. Reproduce Figure 4 setup: re-estimate with an earlier end ----------
estimation.pos.end.fig4 <- idx_to_pos(cal, "2021-04-19")

model <- SSModelDynamicGompertz$new(
  Y = Y, q = q, sea.period = 7,
  start = estimation.pos.start,
  end = estimation.pos.end.fig4,
  calendar = cal
)
res_eval <- estimate(model)

n.ahead <- 14

cat("\n================ FIGURE 4 CASE ================\n")
cat("estimation.end:", tail(res_eval$index, 1), "\n")

# 3. Step through predict_level -----------------------------------------
y.hat.ci <- res_eval$predict_level(n.ahead = n.ahead, sea.on = TRUE,
                                   confidence.level = confidence.level)
cat("y.hat.ci$start:", y.hat.ci$start, "\n")
cat("y.hat.ci positions:", idx_positions(y.hat.ci), "\n")
cat("y.hat.ci$data:\n")
print(y.hat.ci$data)
cat("Any NA/NaN in y.hat.ci$data?:", anyNA(y.hat.ci$data), "\n")

# 4. Step through the eval side, exactly like .plot_holdout_gompertz ----
estimation.end <- tail(res_eval$index, 1)
y.eval.diff <- idx_diff(Y, 1L)

eval_pos <- idx_positions(y.eval.diff)
keep_eval <- eval_pos[eval_pos > estimation.end & eval_pos < estimation.end + n.ahead + 1L]
cat("\nkeep_eval:", keep_eval, "\n")
cat("length(keep_eval):", length(keep_eval), "\n")

y.eval.diff.kept <- y.eval.diff[keep_eval]

common_pos <- intersect(idx_positions(y.eval.diff.kept), idx_positions(y.hat.ci))
cat("common_pos:", common_pos, "\n")
cat("length(common_pos):", length(common_pos), "\n")

d.eval.raw <- data.frame(
  pos = common_pos,
  Actual = idx_values(y.eval.diff.kept[common_pos]),
  Forecast = as.matrix(idx_values(y.hat.ci[common_pos]))[, 1]
)
cat("\nd.eval.raw:\n")
print(d.eval.raw)
cat("anyNA(Actual):", anyNA(d.eval.raw$Actual), "\n")
cat("anyNA(Forecast):", anyNA(d.eval.raw$Forecast), "\n")

d.eval <- na.omit(d.eval.raw)
cat("nrow(d.eval) after na.omit:", nrow(d.eval), "\n")

if (nrow(d.eval) > 0) {
  mape.sea <- mean(100 * (abs(d.eval$Actual - d.eval$Forecast) / d.eval$Actual))
  cat("MAPE:", mape.sea, "\n")
} else {
  cat(">>> d.eval has 0 rows -> mean() of empty vector -> NaN. THIS IS THE BUG SITE. <<<\n")
}

# 5. Also call the real plot_holdout function to confirm the symptom ----
cat("\n================ Calling tsgc::plot_holdout directly ================\n")
p <- tsgc::plot_holdout(
  res = res_eval,
  Y = Y,
  n.ahead = n.forecasts,
  confidence.level = confidence.level,
  series.name = "cases"
)
print(p$labels$subtitle)