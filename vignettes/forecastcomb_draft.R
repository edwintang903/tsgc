# -----------------------------
# 5. Forecast Combination to Predict Hospitalisation from Multiple Lags of Cases
# -----------------------------
# Using the ForecastComb package, combine forecasts for hospitalisations.
Y <- england[, c("cum_cases", "cum_admissions")]
est.start.date <- as.Date("2020-09-01")
est.end.date <- as.Date("2020-10-30")
Y.reinit <- reinitialise_dataframe(Y, est.start.date)

#Plot optimal weights of rolling forecasts
plot_rolling_weights(
  Y.reinit,
  est.start.date,
  est.end.date,
  all_lags = c(2, 5, 7, 9),
  train_days = 20,
  test_days = 60,
  method = comb_BG
)

#Predict future observations with forecast combinations
est.start.date <- as.Date("2020-09-01")
est.end.date <- as.Date("2020-10-30")
idx.est <- (zoo::index(Y.reinit) >= est.start.date) &
  (zoo::index(Y.reinit) <= est.end.date)
y <- Y.reinit[idx.est]

comb_all<-combine_forecasts(
  Y.reinit,
  est.start.date,
  est.end.date,
  all_lags = c(2, 5, 7, 9),
  train_days = 80,
  test_days = 14,
  method = comb_BG
)

#Generate predictions form different lags
future_preds<-matrix(nrow=14,ncol=4)
all_lags<-c(2, 5, 7, 9)
idx.est <- (zoo::index(Y.reinit) >= est.start.date) &
  (zoo::index(Y.reinit) <= est.end.date+80)
y <- Y.reinit[idx.est]
for (i in 1:4){
  j=all_lags[i]
  mod<-SSModelLeadingIndicator(Y=y, n.lag=j)
  resi<-estimate(mod)
  future_preds[,i]<-resi$predict_level(14, sea.on=TRUE)[,1]
}
colnames(future_preds)<-all_lags

#Use weights to combine forecasts
future_combined<-predict(comb_all,future_preds)

#Evaluate the prediction compared to other forecasts
idx.est <- (zoo::index(Y.reinit) >= est.end.date+81) &
  (zoo::index(Y.reinit) <= est.end.date+94)
actual <- diff(Y.reinit)[idx.est,2]

compare<-as.xts(cbind(as.matrix(actual),future_combined, future_preds))
index(compare)<-index(Y.reinit[idx.est,])
colnames(compare)[1]<-"actual"

#mape.comb <- 100*(abs(compare$actual - compare$ForecastTrend)/
#                     compare$Actual) %>% mean %>% round(4)
#mape.sea <- 100*(abs(compare$Actual - compare$Forecast)/compare$Actual) %>%
#  mean %>% round(4)
date_format= "%Y-%m-%d"

df_plot <- as.data.frame(compare)
df_plot$Date <- as.Date(rownames(df_plot), format=date_format)

p1<-ggplot2::ggplot(data = df_plot, aes(x = Date)) +
  ggplot2::geom_line(aes(y = actual, color = "Actual"),lwd = 0.85) +
  ggplot2::geom_line(aes(y = future_combined, color = "Combined Forecast"),lwd = 0.85) +
  geom_line(aes(y = future_preds.2, color = "Lag 2"),lwd = 0.85)+
  ggplot2::geom_line(
    aes(y = future_preds.5, color = "Lag 5"),lwd = 0.85) +
  geom_line(aes(y = future_preds.7, color = "Lag 7"),lwd = 0.85)+
  ggplot2::geom_line(
    aes(y = future_preds.9, color = "Lag 9"),lwd = 0.85) +
  ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045")) +
  ggplot2::geom_ribbon(data = ci_plot, aes(x = Date, ymin = lower, ymax = upper),linetype = 0, linewidth = 0, fill = "#AA2045",
                       alpha = 0.1) +
  labs(x = "Date", y = paste("New",series.name), title = title,
       subtitle = paste("MAPE: ",mape.sea,"%. Trend MAPE: ",
                        mape.trend,"%.",sep="")) +
  theme_economist_white(gray_bg = FALSE, base_size = 14) +
  theme(legend.title = element_blank()) +
  theme(
    text = element_text(size = rel(1)),
    axis.text = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1), margin = margin(r=10)),
    axis.title.x = element_text(size = rel(1), margin = margin(t=10)),
    plot.title = element_text(margin=margin(b=5)),
    plot.subtitle = element_text(
      size = rel(1), hjust=0,  margin = margin(t=3))
  ) +
  scale_linetype_manual(
    values = c("solid", "solid", "solid")) +
  scale_x_date(labels = scales::date_format("%d %b %y")) +
  scale_size_manual(values = c(1, 1.5, 1))