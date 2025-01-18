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
#' @importFrom stats lag
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
#' @param  res.dir File path to save the results to.
#' @param Y Cumulated variable.
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
  y.level.est <- res$data_xts

  # 1. New Cases - Delta Y
  y.hat.diff <- res$predict_level(
    y.cum = y.level.est,
    n.ahead = n.ahead,
    confidence.level= confidence.level,
    sea.on = TRUE,
    return.diff = TRUE
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
  idx.level <- grep("level", colnames(a.t.t))
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



#' @title Returns forecast of number of periods until peak given
#' \code{KFAS::KFS} output.
#'
#' @description Since Harvey and Kattuman (2021) show that \deqn{g_{y,t+\ell|T}
#' = \exp\{\delta_{T|T}+\ell \gamma_{T|T}\}+\gamma_{T|T},} we can compute the
#' \eqn{\ell} for which \eqn{g_{y,t}=0} and then will fall below zero. This
#' \eqn{\ell} is given by
#' \deqn{\ell = \frac{\ln(-\gamma_{T|T})-\delta_{T|T}}{\gamma_{T|T}}.} This is
#' predicated on \eqn{\gamma_{T|T}<0}, else there is super-exponential growth
#' and no peak in sight. Of course, it only makes sense to investigate an
#' upcoming peak for \eqn{g_{y,T|T}>0} (when cases are growing). The estimates
#' of \eqn{\delta_{T|T}} and \eqn{\gamma_{T|T}} are extracted from the
#' \code{KFS} object passed to the function.
#'
#' @param kfs_out The \code{KFAS::KFS} object for which the forecast peak is to
#' be calculated. This would be the \code{output} element of a model estimated
#' in the \code{SSModelDynamicGompertz} or \code{SSModelDynamic}
#'
#' @returns Forecast of number of periods until peak.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' res <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)$estimate()
#'
#' forecast_peak(res$output)
#'
#' @export
forecast_peak <- function(kfs_out) {
  stopifnot(class(kfs_out) == "KFS")
  n <- attr(kfs_out$model, "n")
  delta.t.t <- kfs_out$att[n, 'level']
  gamma.t.t <- kfs_out$att[n, 'slope']
  if (gamma.t.t > 0) {
    return(Inf)
  } else {
    return(forecast.peak(delta = delta.t.t, gamma = gamma.t.t)
    )
  }
}


#' @title Returns forecast of number of periods until peak given estimated
#' state variables \eqn{\delta} and \eqn{\gamma}.
#'
#' @description Since Harvey and Kattuman (2021) show that
#' \deqn{g_{y,t+\ell|T} = \exp\{\delta_{T|T}+\ell \gamma_{T|T}\}+\gamma_{T|T},}
#' we can compute the \eqn{\ell} for which \eqn{g_{y,t}=0} and then will fall
#' below zero. This \eqn{\ell} is given by
#' \deqn{\ell = \frac{\ln(-\gamma_{T|T})-\delta_{T|T}}{\gamma_{T|T}}.} This is
#' predicated on \eqn{\gamma_{T|T}<0}, else there is super-exponential growth an
#' no peak in sight. Of course, it only makes sense to investigate an upcoming
#' peak for \eqn{g_{y,T|T}>0} (when cases are growing).
#'
#' @param delta The estimate of \eqn{\delta}, the level of \eqn{\ln g}.
#' @param gamma The estimate of \eqn{\gamma}, the slope of \eqn{\ln g}.
#'
#' @examples
#' # Forecasts the peak of an epidemic with gamma < 0 so that a peak is in
#' # sight.
#' forecast.peak(-2.87,-0.045)
#'
#' # Does not return a result (returns an error as gamma > 0)
#' try(forecast.peak(-2.87,0.045), silent=TRUE)
#'
#' @returns Forecast of number of periods until peak.
#'
#' @export
forecast.peak <- function(delta, gamma) {
  # if numerator positive then get negative forecast for ell^*
  # stopifnot(log(-gamma) - delta < 0)
  # If gamma positive then no saturation level.
  stopifnot(gamma < 0)
  return(
    as.numeric(
      (log(-gamma) - delta) / gamma
    ))
}

#' @title Compute Mean Absolute Percentage Error (MAPE) for Forecasts Against
#' a Holdout Sample
#'
#' @description This is a helper function that calculates the Mean Absolute
#' Percentage Error (MAPE) of a forecast generated by time series growth curve
#' (tsgc) models. It compares the forecasted values to a holdout sample,
#' providing a measure of forecast accuracy. If the vanilla growth curve model
#' is used, ensure that the date format is "%Y-%m-%d".
#'
#' @param object A `filterResults` or `filterResultsLI` object, obtained from
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
mapes<-function(object,...){
  object$mapes(...)
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
#' est.end.date=estimation.date.end,n.ahead=7,lags=1:9,totaldays=3, 
#' vanilla=TRUE,freq=2,LeadIndCol=1)
#'
#' @export
cross_val<-function(y,est.end.date,n.ahead,all_lags,totaldays=1,freq=1, vanilla=TRUE,
                    LeadIndCol=1){
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
      res <- model_q$estimate()
      results[1,k+1]=round(mapes(res,n.ahead,Z)$sea,2)
    }
    for (i in all_lags){
      out<-SSModelLeadingIndicator(Y=y[index(y) <= est.end.date+(k-1)*freq],n.lag = i)
      res<-out$estimate()
      results[vanilla+i,k+1]<-round(mapes(res,n.ahead,y)$sea,2)
    }
    results[length(all_lags)+vanilla+1,k+1]=allall_lags[which.min(results[,k+1])]
  }
  alldates<-as.character(est.end.date+c(0:(k-1))*freq)
  colnames(results)<-c("Lag",alldates)
  return(results)
}


#' @title Optimal Forecast Combination for Dynamic Gompertz and Leading Indicator Models
#'
#' @description This function computes the optimal combination of forecasts from 
#' the Dynamic Gompertz model and the Leading Indicator model over specified lags. 
#' It provides flexibility to train and test the combination model, optionally 
#' displaying the dynamic evolution of weights over time.
#'
#' @param Y An `xts` object with a date index representing the cumulative data series. 
#' The object should contain two columns: the leading indicator and the target variable. 
#' The specific column for the leading indicator can be specified using the \code{LeadIndCol} parameter.
#' 
#' @param est.start.date A `Date` object indicating the start date of the base estimation dataset.
#' Forecasts generated by the model will include data starting from this date.
#' 
#' @param est.end.date A `Date` object indicating the end date of the base estimation dataset.
#' The model will use data within the range defined by \code{est.start.date} and \code{est.end.date}
#' for initial training and estimation. Predictions will then extend beyond \code{est.end.date}
#' over a period determined by the sum of \code{train_days} and \code{test_days}. Specifically:
#' - The period \code{est.end.date + 1 : train_days} is used for training the forecast combination model.
#' - The period \code{est.end.date + train_days + 1 : test_days} is used for testing and evaluation.
#' 
#' @param all_lags A numeric vector specifying the lags to be used in the 
#' Leading Indicator model. A 0 value indicates the Dynamic Gompertz Model is used.
#' 
#' @param train_days An integer specifying the number of days used to train the 
#' forecast combination model.
#' 
#' @param test_days An integer specifying the number of days used to test the 
#' forecast combination model.
#' 
#' @param method A function or method from the \code{ForecastComb} package to 
#' compute the optimal 
#' combination of forecasts. Example: \code{comb_OLS}.
#' 
#' @param LeadIndCol An integer indicating the column number in \code{Y} that 
#' contains the leading indicator. Defaults to 1.
#' 
#' @param rolling A logical value. If TRUE, computes the forecast combination 
#' dynamically, showing the evolution of weights across \code{test_days}. 
#' Defaults to FALSE.
#' 
#' @importFrom ggplot2 ggplot geom_line geom_point labs theme_minimal
#' @importFrom tidyr pivot_longer
#' @import ForecastComb
#'
#' @return If \code{rolling = TRUE}, returns a list with two elements:
#' - A forecast combination object with details of the computed weights.
#' - A `ggplot` object showing the evolution of weights over time.
#' 
#' If \code{rolling = FALSE}, returns the results of the forecast combination 
#' method applied to the data.
#'
#' @examples
#' library(tsgc)
#' library(ForecastComb)
#' library(tidyr)
#' Y <- england[, c("cum_cases", "cum_admissions")]
#' est.start.date <- as.Date("2020-09-01")
#' est.end.date <- as.Date("2020-10-30")
#' Y.reinit <- reinitialise_dataframe(Y, est.start.date)
#' 
#' combine_forecasts(
#'   Y.reinit, 
#'   est.start.date, 
#'   est.end.date,
#'   all_lags = c(2, 5, 7, 9), 
#'   train_days = 20, 
#'   test_days = 60,
#'   method = comb_BG, 
#'   rolling = TRUE
#' )
#'
#' @export
combine_forecasts=function(Y,est.start.date,est.end.date,
                           all_lags,train_days,test_days,method,
                           LeadIndCol=1,rolling=FALSE){  
  
  x <- Value <- Variable <- NULL
  
  if (!inherits(Y, "xts")) stop("Y must be an xts object.")
  if (!is.numeric(all_lags) || any(all_lags < 0)) stop("all_lags must be a numeric vector with non-negative values.")
  if (!is.function(method)) stop("method must be a valid function from the ForecastComb package.")
  if (!is.logical(rolling)) stop("rolling must be a logical value.")
  if (!is.numeric(LeadIndCol) || LeadIndCol <= 0) stop("LeadIndCol must be a positive integer.")
  
  Y_full<-add_daily_ldl(Y,LeadIndCol=LeadIndCol)
  lag.length=length(all_lags)
  
  #calculating all predictions
  total_days=train_days+test_days
  
  result <- matrix(0, nrow = lag.length, ncol = total_days)
  actual <- matrix(0, nrow = 1, ncol = total_days)
  for (k in 1:lag.length){
    for (m in c(0:(total_days-1))){
      idx<-(index(Y) <= est.end.date+m) & (index(Y) >= est.start.date)
      if (all_lags[k]==0){
        out<-SSModelDynamicGompertz(Y=Y[idx,-LeadIndCol])
        res<-out$estimate()
        d<-res$predict_level(y.cum=Y[idx,-LeadIndCol],n.ahead=1,confidence.level = 0.68, 
                             return.diff = TRUE)
        result[k,1+m]<-as.matrix(d)[1,1]
      } else{
        out<-SSModelLeadingIndicator(Y=Y[idx], n.lag = all_lags[k])
        res<-out$estimate()
        result[k,1+m]<-res$predict_level(n.ahead=1)$seasonal[,"forc"]
      }
      actual[1+m]<-Y_full[est.end.date+m+1, "newAdmit"]
    }
  }
  train_o<-actual[1:train_days]
  train_p<-t(result[,1:train_days])
  test_o<-actual[(train_days+1):total_days]
  test_p<-t(result[,(train_days+1):total_days])
  data<-foreccomb(train_o, train_p, test_o, test_p)
  
  if (rolling){
    combin<-rolling_combine(data, deparse(substitute(method)))
    df<-as.data.frame(combin$Weights)
    colnames(df)<-all_lags
    df$x<-est.end.date+index(df)+train_days
    
    # Reshape the data from wide to long format
    df_long <- df %>%
      pivot_longer(names_to= "Variable", values_to = "Value", cols=-x)
    
    # Plot the data
    p1<-ggplot(df_long, aes(x = x, y = Value, color = Variable, group = Variable)) +
      geom_line(linewidth = 1) +
      geom_point(size= 2) +
      labs(title = "Evolution of weights",
           x = "Observations",
           y = "Values") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(list(combin, p1))
  } else {
    return(method(data))
  }
}




