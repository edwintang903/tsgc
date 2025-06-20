# Created by: Craig Thamotheram
# Created on: 19/02/2022

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' @title Compute log growth rate of cumulated dataset
#
#' @description Helper method to compute the log growth rates of cumulated
#' variables. It will compute the log cumulative growth rate for each column in
#' the data frame.
#'
#' @param dt Cumulated data series.
#' @returns A data frame of log growth rates of the cumulated variable which has
#' been inputted via the parameter \code{dt}.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' df2ldl(gauteng)
#'
#'
#' @export
df2ldl <- function(dt) {
  dt.ldl <- log(diff(dt) / stats::lag(dt))
  return(dt.ldl)
}

#' @title Subsetting xts objects given start dates and end dates
#
#' @description Helper method to subset a data frame for a specified time frame
#'
#' @param df xts object
#' @param start.date Start date of time frame, Date object.
#' @param end.date End date of time frame, Date object.
#' @returns A subsetted data frame
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#'
#'
#' @export
get_timeframe<-function(df, start.date, end.date=NULL){
  if (is.null(end.date)){
    idx.est1 <- (zoo::index(df) >= start.date)
  } else {
    idx.est1 <- (zoo::index(df) >= start.date) & (zoo::index(df) <= end.date)
  }
  return(df[idx.est1,])
}


#' @title Compute successive increments and log growth rate of 2-variable
#' cumulated dataset
#
#' @description Helper method to compute the successive increments and log
#' growth rates of cumulated variables. It will compute the successive
#' increments and log cumulative growth rate for each column in the
#' 2-column data frame, which will then be used to predict or estimate with the
#' leading indicator model.
#'
#' @param data Cumulated data series in xts format with date index and columns:
#' leading indicator and target variable. Can specify which column is leading
#' indicator by \code{LeadIndCol} parameter.
#' @param LeadIndCol Column number of \code{data} that contains the leading
#' indicator
#' @returns A data frame with original cumulative variable, successive increments
#' and log growth rates.
#'
#' @examples
#' library(tsgc)
#' data(england,package="tsgc")
#' add_daily_ldl(england[,c("cum_cases","cum_admissions")],LeadIndCol=1)
#'
#' @importFrom xts lag.xts
#'
#' @export
add_daily_ldl <- function(data, LeadIndCol=1){
  if (LeadIndCol==1){
    names(data)<-c("cCases", "cAdmit")
  } else {
    names(data)<-c("cAdmit", "cCases")
    data<-data[,c(2,1)]
  }
  data$newCases = diff(data$cCases)
  data$newAdmit = diff(data$cAdmit)
  data$LDLcases = log(as.vector(data$newCases)/lag(as.vector(data$cCases)))
  data$LDLhosp = log(as.vector(data$newAdmit)/lag(as.vector(data$cAdmit)))
  return(data)
}


#' @title Reinitialise a data frame by subtracting the `reinit.date` row from
#' all columns
#'
#' @param dt Cumulated data series.
#' @param reinit.date Reinitialisation date. E.g. \samp{'2021-05-12'}.
#'
#' @returns The reinitialised data frame
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' reinitialise_dataframe(gauteng,as.Date("2021-01-01"))
#'
#' @export
reinitialise_dataframe <- function(dt, reinit.date) {
  # Take cumulative dataframe and reinit from reinit.date as first date of data
  # 1. Get data frame including date before reinit.date
  dt <- dt[index(dt) >= as.Date(reinit.date) - 1,]

  # 2. Substract away the t-1 date data
  dt <- sweep(dt, 2, dt[1,])

  # 3. Keep only data from t onwards.
  dt <- dt[index(dt) >= as.Date(reinit.date),]
  return(dt)
}


#' @title Return index and value of maximum
#' @description Something similar to Python's argmax.
#' @param x Object to have its maximum found
#' @param decreasing Logical value indicating whether \code{x} should be
#' ordered in decreasing order. Default is \code{TRUE}. Setting this to
#' \code{FALSE} would find the minimum.
#' @returns The maximum value and its index.
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' argmax(gauteng)
#' @export
argmax <- function(x, decreasing=TRUE) {
  return(x[order(x, decreasing = decreasing)[1]])
}


#' @title Write a selection of relevant results to disc
#'
#' @description Function writes the following results to csv files which get
#' saved in the location specified in \code{res.dir}: forecast new cases or
#' incidence variable, \eqn{y}; the filtered level and slope of \eqn{\ln g},
#' \eqn{\delta} and \eqn{\gamma}; filtered estimates of \eqn{g_y} and the
#' confidence intervals for these estimates.
#'
#' @param res Results object estimated using the \samp{estimate()} method.
#' @param res.dir File path to save the results to.
#' @param n.ahead Number of periods ahead to forecast.
#' @param confidence.level Confidence level to use for the confidence interval
#' on the forecasts \eqn{\ln(g_t)}.
#' 
#' @importFrom utils write.csv
#' @importFrom stats qnorm
#'
#' @returns A number of csv files saved in the directory specified in
#' \code{res.dir}.
#' @examples
#' # Not run as do not wish to save to local disc when compiling documentation.
#' # Below will run if copied and pasted into console.
#' library(tsgc)
#' library(here)
#'
#' res.dir <- tempdir()
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#' res <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)$estimate()
#'
#' tsgc::write_results(
#' res=res, res.dir = res.dir, Y = gauteng[idx.est], n.ahead = 14,
#' confidence.level = 0.68
#' )
#'
#' @export
write_results <- function(res, res.dir, n.ahead, confidence.level=0.68) {
  # 1. New Cases - Delta Y
    y.hat.diff <- res$predict_level(
      n.ahead = n.ahead,
      confidence.level= confidence.level,
      sea.on = TRUE
    )
  
  write.csv(
    y.hat.diff,
    row.names = index(y.hat.diff),
    file = file.path(res.dir, "y-forecast.csv")
  )
  
  # 2. Filtered slope / level
  y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
  filtered.level <- y.hat.all$level.t.t
  filtered.slope <- y.hat.all$slope.t.t
  a.t.t <- y.hat.all$a.t.t
  P.t.t <- y.hat.all$P.t.t
  idx.slope <- grep("slope", colnames(a.t.t))
  idx.level <- grep("level", colnames(a.t.t))[1]
  gamma.std.err <- sqrt(P.t.t[idx.slope, idx.slope,])
  delta.std.err <- sqrt(P.t.t[idx.level, idx.level,])
  gamma <- cbind(filtered.slope, gamma.std.err)
  delta <- cbind(filtered.level, delta.std.err)
  colnames(gamma) <- c("gamma", "std.err")
  colnames(delta) <- c("delta", "std.err")
  write.csv(
    gamma,
    row.names = index(filtered.slope),
    file = file.path(res.dir, "gamma_filtered.csv")
  )
  write.csv(
    delta,
    row.names = index(filtered.level),
    file = file.path(res.dir, "delta_filtered.csv")
  )

  # 3. Filtered growth rate of new cases (g_{y}) - CI from standard error on
  # slope component of
  # state covariance matrix.
  g.y.t.t <- exp(filtered.level) + filtered.slope
  ci <- qnorm((1 - confidence.level) / 2) * gamma.std.err %o% c(1, -1)
  ci_bounds <- as.vector(g.y.t.t) + ci
  gy.ci <- xts(cbind(fit = g.y.t.t, prediction = ci_bounds),
               order.by = index(filtered.level))
  names(gy.ci)[2:3] <- list('lower', 'upper')

  write.csv(
    gy.ci,
    row.names = index(g.y.t.t),
    file = file.path(res.dir, "g_y_filtered.csv")
  )

}

#' @title Calculate reproduction number estimates and credible intervals
#'
#' @description Following Harvey and Kattuman (2021)
#'
#' @param res A `filterResults` object, obtained from \code{estimate()} method.
#' @param gen_int Generation interval in days
#' @param ndays Number of days to plot, counting from the end of estimation timeframe.
#' @param show_plot A logical value indicating whether ti show the plot of R0
#' @param title Title for the reproduction number plot. 
#'
#' @returns Graph of estimated Rt and forecast intervals.
#' 
#' @importFrom timetk tk_tbl
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#'
#' @export
estimate_r0<-function(res,gen_int, ndays=7, show_plot=FALSE, 
                      title="Reproduction number"){
  r.t <- tail(exp(res$get_gy_ci() * gen_int), ndays) %>% tk_tbl
  names(r.t) <- c("Date", "Rt", "lower", "upper")
  
  if (show_plot){
    res.rt <- ggplot(r.t, aes(x = Date)) +
      ylim(0, 1.4) +
      geom_line(aes(y = Rt, color = "Rt")) +
      geom_point(aes(y = Rt), color = "red", size = 3) +
      geom_segment(aes(xend = Date, yend = lower, y = Rt), color = "blue") +
      geom_segment(aes(xend = Date, yend = upper, y = Rt), color = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = "68%  Interval"), alpha = 0.2) +
      geom_hline(yintercept = 1, linetype = "solid", linewidth = 1.5, color = "black") +
      scale_x_date(date_breaks = "1 day") +
      labs(title = title)+
      theme_light(base_size = 12) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(face = "bold")
      )
    return(res.rt)
  } else {
    return(r.t)
  }
}

#' @title Compute Mean Absolute Percentage Error (MAPE) for Forecasts Against
#' a Holdout Sample
#'
#' @description This is a helper function that calculates the Mean Absolute
#' Percentage Error (MAPE) of a forecast generated by time series growth curve
#' (tsgc) models. It compares the forecasted values to a holdout sample,
#' providing a measure of forecast accuracy. 
#'
#' @param res A `filterResults` or `filterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param n.ahead Integer specifying the number of days to forecast ahead.
#' @param Y An xts object containing the original cumulative dataset.
#'
#' @returns A list containing the MAPE values for both the trend forecast and
#' the forecast that includes the seasonal component.
#' 
#' @examples
#' library(tsgc)
#' 
#' #Setup
#' date_format="%Y-%m-%d"
#' estimation.date.start=as.Date("2021-04-30")
#' estimation.date.end=as.Date("2021-07-24")
#' n.ahead=7
#' Y=gauteng
#' idx.est =(zoo::index(Y) >= estimation.date.start) & (zoo::index(Y) <= estimation.date.end)
#' y = Y[idx.est]
#' 
#' #Estimate the model
#' model_q <- SSModelDynamicGompertz$new(Y = y)
#' res <- estimate(model_q)
#' 
#' #Return MAPE of forecast
#' mapes(res,n.ahead=n.ahead,Y)
#' 
#' 
#' @export
mapes<-function(res,n.ahead,Y){
  res$mapes(n.ahead,Y)
}


#' @title Walk-Forward Validation for Lag Comparison Using Mean Absolute
#' Percentage Error (MAPE)
#'
#' @description This function performs a walk-forward validation to compare
#' forecasting performance across different lag values. It returns a table of
#' MAPE values for forecasts \code{n.ahead} days ahead, using models estimated
#' with varying all_lags over a series of specified end dates.
#'
#' @param y An xts object representing the cumulative data series with a date
#' index. The object should include columns for both the leading indicator and
#' the target variable. The specific column for the leading indicator can be
#' designated using the \code{LeadIndCol} parameter.
#' @param est.end.date The initial estimation end date for model fitting.
#' Starting from this date, the function re-estimates the model and evaluates
#' the performance for each lag in \code{all_lags} every \code{freq} days, over a
#' period of \code{totaldays} days.
#' @param n.ahead Integer specifying the number of days to forecast ahead for
#' MAPE evaluation.
#' @param all_lags Positive integer-valued array specifying which lags we are using
#' in leading indicator models for comparison.
#' @param totaldays Integer indicating the total number of days for which
#' walk-forward validation results will be reported.
#' @param freq Integer specifying the frequency, in days, at which the model
#' is re-estimated and evaluated during the walk-forward validation.
#' @param vanilla Logical. If \code{TRUE}, the function compares the vanilla
#' growth curve model to the leading indicator models with different all_lags. The
#' results for the vanilla model are presented in the row where Lag=0.
#' Note: If using the vanilla model, the date format must be "%Y-%m-%d".
#' @param LeadIndCol (Only required for leading indicator models) Integer
#' representing the column number in \code{y} that contains the leading
#' indicator.
#' @param criterion A string object indicating how to compare between different 
#' models. Available choices are "mape" (by default), "mae" and "rmse". 
#' 
#' @importFrom zoo index
#'
#' @returns A table summarizing the MAPE scores for each lag across the
#' specified dates.
#'
#' @examples
#' library(tsgc)
#' 
#' #Lay out the estimation settings
#' Y = england[,1:2] 
#' estimation.date.start = as.Date("2021-04-30")
#' estimation.date.end = as.Date("2021-07-24")
#' 
#' #Output cross validation result
#' cross_val(y=Y[index(Y)>=estimation.date.start],
#' est.end.date=estimation.date.end,n.ahead=7,all_lags=1:9,totaldays=3, 
#' vanilla=TRUE,freq=2,LeadIndCol=1, criterion="mae")
#'
#' @export
cross_val<-function(y,est.end.date,n.ahead,all_lags,totaldays=1,freq=1, vanilla=TRUE,
                    LeadIndCol=1, criterion="mape"){
  if (vanilla){
    allall_lags<-c(0,all_lags)
  }
  else{
    allall_lags<-all_lags
  }
  results <- data.frame(
    Lag = c(allall_lags,"Min MAPE at")
  )
  for (k in 1:totaldays){
    if (vanilla){
      date_format="%Y-%m-%d"
      Z = y[,-LeadIndCol]
      model_q <- SSModelDynamicGompertz$new(Y = Z[index(Z) <= est.end.date+(k-1)*freq])
      res <- estimate(model_q)
      results[1,k+1]=round(mapes(res,n.ahead,Z)[[criterion]],2)
    }
    for (i in all_lags){
      out<-SSModelLeadingIndicator(Y=y[index(y) <= est.end.date+(k-1)*freq],n.lag = i)
      res<-estimate(out)
      results[vanilla+i,k+1]<-round(mapes(res,n.ahead,y)[[criterion]],2)
    }
    results[length(all_lags)+vanilla+1,k+1]=allall_lags[which.min(results[,k+1])]
  }
  alldates<-as.character(est.end.date+c(0:(k-1))*freq)
  colnames(results)<-c("Lag",alldates)
  return(results)
}