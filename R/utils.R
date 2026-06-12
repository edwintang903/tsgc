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
#' # Select data between 2020-04-01 and 2020-08-01, inclusive
#' get_timeframe(gauteng, as.Date("2020-04-01"), as.Date("2020-08-01"))
#' 
#' # Select all data after 2020-04-01, inclusive
#' get_timeframe(gauteng, as.Date("2020-04-01"))
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
#' @param data Cumulated data series in xts format with date index and 2 columns:
#' leading indicator and target variable. Can specify which column is leading
#' indicator by \code{LeadIndCol} parameter.
#' @param LeadIndCol Column number of \code{data} that contains the leading
#' indicator. An integer that can only take values 1 (by default) or 2.
#' @returns A data frame with original cumulative variable, successive increments
#' and log growth rates.
#'
#' @examples
#' library(tsgc)
#' data(england,package="tsgc")
#' add_daily_ldl(england[,c("cum_cases","cum_admissions")],LeadIndCol=1)
#'
#' @importFrom xts is.xts
#'
#' @export
add_daily_ldl <- function(data, LeadIndCol=1){
  if (!is.xts(data)){
    stop("data is not an xts object.")
  }
  if (LeadIndCol==1){
    names(data)<-c("cLead", "cTarg")
  } else if (LeadIndCol==2) {
    names(data)<-c("cTarg", "cLead")
    data<-data[,c(2,1)]
  } else {
    stop("LeadIndCol must be an integer, either 1 or 2.")
  }
  ldl <- df2ldl(data)
  
  data$newLead = diff(data$cLead)
  data$newTarg = diff(data$cTarg)

  data$LDLlead = ldl$cLead
  data$LDLtarg = ldl$cTarg
  return(data)
}


#' @title Reinitialise a data frame by subtracting the `reinit.date` row from
#' all columns
#'
#' @param dt Cumulated data series, belonging to the xts class. Must only contain 
#' 1 data column in addition to a date index.
#' @param reinit.date Reinitialisation date, belonging to Date, yearmon or yearqtr classes. E.g. \samp{as.Date('2021-05-12')}.
#'
#' @returns The reinitialised data frame
#' @importFrom zoo index
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' reinitialise_dataframe(gauteng,as.Date("2021-01-01"))
#'
#' @export
reinitialise_dataframe <- function(dt, reinit.date) {
  if (dim(dt)[2]!=1){
    stop("dt must only contain 1 data column in addition to a date column.")
  }
  # Take cumulative dataframe and reinit from reinit.date as first date of data
  # 1. Get data frame including date before reinit.date
  first_ind<-which.max(index(dt) == reinit.date)
  dt <- dt[(first_ind-1):length(dt),]

  # 2. Subtract away the t-1 date data
  dt <- sweep(dt, 2, dt[1,])

  # 3. Keep only data from t onwards.
  dt <- dt[-1,]
  return(dt)
}


#' @title Return index and value of maximum
#' @description Similar to Python's argmax function.
#' @param x Object to have its maximum found, usually an xts object.
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


#' @title Write a selection of relevant results to disk
#'
#' @description Function writes the following results to csv files which get
#' saved in the location specified in \code{res.dir}: forecast new cases or
#' incidence variable, \eqn{y}; the filtered level and slope of \eqn{\ln g},
#' \eqn{\delta} and \eqn{\gamma}; filtered estimates of \eqn{g_y} and the
#' confidence intervals for these estimates.
#'
#' @param res Results object of class \code{FilterResults} or \code{FilterResultsLI}, 
#' obtained from \samp{estimate()} method.
#' @param res.dir File path to save the results to. A character string.
#' @param n.ahead Number of periods ahead to forecast. A positive integer.
#' @param prefix The prefix to be added to the file names generated. A character string. 
#' @param confidence.level Confidence level to use for the confidence interval
#' on the forecasts \eqn{\ln(g_t)}.
#' 
#' @importFrom utils write.csv
#' @importFrom stats qnorm
#'
#' @returns A number of csv files saved in the directory specified in
#' \code{res.dir}.
#' @examples
#' # Not run as do not wish to save to local disk when compiling documentation.
#' # Below will run if copied and pasted into console.
#' library(tsgc)
#' library(here)
#'
#' res.dir <- tempdir()
#' data(gauteng,package="tsgc")
#' res <- estimate(SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, 
#' end.date=as.Date("2020-07-06")))
#'
#' tsgc::write_results(
#' res=res, res.dir = res.dir, prefix="dyn_gompertz",n.ahead = 14,
#' confidence.level = 0.68)
#'
#' @export
write_results <- function(res, res.dir, n.ahead, prefix="", confidence.level=0.68) {
  if (class(res)!="FilterResults" && class(res)!="FilterResultsLI"){
    stop("res must be a FilterResults or FilterResultsLI object.")
  }
  # 1. New Cases - Delta Y
    y.hat.diff <- res$predict_level(
      n.ahead = n.ahead,
      confidence.level= confidence.level,
      sea.on = TRUE)
  
  write.csv(
    y.hat.diff,
    row.names = index(y.hat.diff),
    file = file.path(res.dir, paste(prefix,"cases_fcst.csv", sep="")))
  
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
    file = file.path(res.dir, paste(prefix, "trend_slope_filt.csv", sep=""))
  )
  write.csv(
    delta,
    row.names = index(filtered.level),
    file = file.path(res.dir, paste(prefix, "log_gr_level_filt.csv", sep=""))
  )

  # 3. Filtered growth rate of new cases (g_{y}) - CI from standard error on
  # slope component of state covariance matrix.
  g.y.t.t <- exp(filtered.level) + filtered.slope
  ci <- qnorm((1 - confidence.level) / 2) * gamma.std.err %o% c(1, -1)
  ci_bounds <- as.vector(g.y.t.t) + ci
  gy.ci <- xts(cbind(fit = g.y.t.t, prediction = ci_bounds),
               order.by = index(filtered.level))
  colnames(gy.ci)[2:3] <- c('lower', 'upper')

  write.csv(
    gy.ci,
    row.names = index(g.y.t.t),
    file = file.path(res.dir, paste(prefix, "cases_gr.csv", sep="")))
  
  message("Saved results for: ", substitute(res))
  
}

#' @title Calculate reproduction number estimates and credible intervals
#'
#' @description Reproduction number is estimated based on the method described in Harvey and Kattuman (2021).
#'
#' @param res A `FilterResults` object, obtained from \code{estimate()} method.
#' @param gen_int Generation interval in days
#' @param ndays Number of days to plot, counting from the end of estimation timeframe.
#' @param show_plot A logical value indicating whether to show the plot of R0
#' @param title Title for the reproduction number plot. 
#'
#' @returns Graph of estimated Rt and forecast intervals.
#' 
#' @importFrom timetk tk_tbl
#' @import ggplot2
#' 
#' @references Harvey A, Kattuman P (2021). “A farewell to R: time-series models 
#' for tracking and forecasting epidemics.” Journal of the Royal Society Interface, 18. URL http://doi.org/10.1098/rsif.2021.0179.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' cumulative_cases <- gauteng[, 1] 
#' 
#' # Estimate model
#' model_q <- SSModelDynamicGompertz$new(Y = cumulative_cases, q = 0.005,
#'                                       start.date=as.Date("2021-02-01"), 
#'                                       end.date=as.Date("2021-04-19"))
#' res_q <- estimate(model_q)
#' summary(res_q)
#' 
#' # Calculate reproduction number estimates and credible intervals.
#' gen_int <- 4  # Generation interval in days
#' ndays<-7 #Number of days to plot
#' r.t <- estimate_r0(res_q, gen_int, ndays)
#' r.t 
#' 
#'  # Plot reproduction numbers.
#' estimate_r0(res_q, gen_int, ndays, show_plot = TRUE, title="Gauteng Reproduction numbers")
#'
#' @export
estimate_r0<-function(res, gen_int, ndays=7, show_plot=FALSE, 
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
#' @description This is a helper function that calculates five error metrics of 
#' a forecast generated by time series growth curve
#' (tsgc) models. It compares the forecasted values to a holdout sample,
#' providing a measure of forecast accuracy. 
#'
#' @param res A `FilterResults` or `FilterResultsLI` object, obtained from
#' \code{estimate()} method.
#' @param n.ahead Integer specifying the number of days to forecast ahead.
#' @param Y An xts object containing the original cumulative dataset.
#'
#' @returns A list containing five error metrics for the forecast, with element names
#' \itemize{
#' \item mape: mean absolute percentage error 
#' \item smape: symmetric mean absolute percentage error (between 0 to 100) 
#' \item mae: mean absolute error 
#' \item rmse: root mean squared error 
#' \item coverage: Percentage of holdout sample data points that lie inside the 
#' confidence interval for predictions}
#' 
#' @examples
#' library(tsgc)
#' 
#' #Estimate the model
#' model_q <- SSModelDynamicGompertz$new(Y = gauteng, start.date=as.Date("2021-04-30"), end.date=as.Date("2021-07-24"))
#' res <- estimate(model_q)
#' 
#' #Return MAPE of forecast
#' mapes(res,n.ahead=7,gauteng)
#' 
#' 
#' @export
mapes<-function(res,n.ahead,Y){
  res$mapes(n.ahead,Y)
}


#' @title Walk-Forward Validation for Model Comparison Using Mean Absolute
#' Percentage Error (MAPE)
#'
#' @description This function performs a walk-forward validation to compare
#' forecasting performance across different models specified by the user. 
#' It returns a data frame of a user-specified error metric (e.g. MAPE, MAE) 
#' for forecasts \code{n.ahead} days ahead, using the given models with varying
#' end dates.
#'
#' @param Y An xts object representing the cumulative data series with a date
#' index. If a Leading Indicator model is compared, Y should include columns for both the leading indicator and
#' the target variable. The specific column for the leading indicator can be
#' designated using the \code{LeadIndCol} parameter.
#' @param model_list A list containing \code{SSModelDynamicGompertz} or \code{SSModelLeadingIndicator}
#' objects, to be compared in a cross validation procedure.
#' @param est.end.date The initial estimation end date for model fitting.
#' Starting from this date, the function re-estimates the model and evaluates
#' the performance for each lag in \code{all_lags} every \code{gap} days, over a
#' period of \code{n.estimate} days.
#' @param n.ahead Integer specifying the number of days to forecast ahead for
#' MAPE evaluation.
#' @param n.estimate Integer indicating the total number of days for which
#' walk-forward validation results will be reported.
#' @param gap Integer specifying the time gap between two successive validations, where the model
#' is re-estimated and evaluated during the walk-forward validation.
#' @param xpred_lead.full (Only for required for leading indicator models) 
#' An xts object containing the values of exogenous variables for 
#' the leading indicator over the estimation and prediction time frame.
#' @param xpred_targ.full An xts object containing the values of exogenous variables for 
#' the target variable over the estimation and prediction time frame.
#' @param LeadIndCol (Only required for leading indicator models) Integer
#' representing the column number in \code{y} that contains the leading
#' indicator.
#' @param criterion A string object indicating how to compare between different 
#' models. Available choices are "mape" (by default), "smape", "mae" and "rmse". 
#' 
#' @importFrom zoo index
#' @importFrom magrittr and
#'
#' @returns A table summarizing the chosen error metric for each model in
#' \code{model_list} across the specified dates.
#'
#' @examples
#' library(tsgc)
#' library(KFAS)
#' 
#' #Lay out the estimation settings
#' est.start <- as.Date("2020-02-25")
#' est.end <- as.Date("2020-04-01")
#' Yuk <- tsgc::ukitaly[, "UK"]
#' 
#' # Example 1: ukitaly dataset 
#' # Create a list to store different models
#' cv_models<-list()
#' 
#' # Model 1: Vanilla Gompertz
#' cv_models[["Vanilla_q"]]<-SSModelDynamicGompertz(Y=Yuk, q=0.005, start.date = est.start, end.date = est.end)
#' 
#' # Model 2: Vanilla Gompertz with AR1
#' cv_models[["Vanilla_ar1"]]<-SSModelDynamicGompertz(Y=Yuk, start.date = est.start, end.date = est.end, ar1=TRUE)
#' 
#' # Model 3-6: Leading Indicator with different lags from 7, 10, 14 or 18
#' for (i in c(7,10,14,18)){
#'   cv_models[[paste0("Lag", i)]]<-SSModelLeadingIndicator(Y=ukitaly, start.date = est.start, end.date = est.end, n.lag=i)}
#' 
#' # Display cross-validation analysis
#' cross_val(Y=ukitaly, model_list=cv_models, est.end.date = est.end, n.estimate=5, gap=2)
#'
#' # Example 2: England hospitalizations (with xpred)
#' eng <- tsgc::england[, 1:2]
#' est.start.eng <- as.Date("2021-04-30")
#' est.end.eng   <- as.Date("2021-07-24")
#' 
#' #Cross-validation example 
#' cv_models=list()
#' # Model 1: Vanilla Gompertz
#' cv_models[["Vanilla_q"]]<-SSModelDynamicGompertz(Y=eng[,-1], q=0.005, start.date = est.start, end.date = est.end)
#' 
#' # Model 2: Vanilla Gompertz with xpred
#' cv_models[["Vanilla_xpred"]]<-SSModelDynamicGompertz(Y=eng[,-1], start.date = est.start.eng, end.date = est.end.eng, xpred=england_weather_2021)
#' 
#' # Model 3-6: Leading Indicator with different lags or with xpred
#' for (i in c(3,4)){
#' cv_models[[paste0("Lag", i)]]<-SSModelLeadingIndicator(Y=eng, 
#' start.date = est.start.eng, end.date = est.end.eng, n.lag=i)
#' 
#' cv_models[[paste0("Lag", i,"_xpred")]]<-SSModelLeadingIndicator(Y=eng, 
#' start.date = est.start.eng, end.date = est.end.eng, 
#' n.lag=i, xpred_lead = england_weather_2021, xpred_targ=england_weather_2021)}
#' 
#' # Display cross-validation analysis
#' cross_val(Y=eng, model_list=cv_models, est.end.date = est.end.eng, n.estimate=5, 
#' gap=2, xpred_targ.full = england_weather_2021, xpred_lead.full = england_weather_2021)
#'
#' @export
cross_val<-function(Y, model_list, est.end.date, n.ahead=7, n.estimate=1, gap=1,
                    xpred_targ.full=NULL, xpred_lead.full=NULL, LeadIndCol=1, criterion="mape"){
  if (!is.xts(Y)){
    stop("Y must be an xts object.")
  }
  if (dim(Y)[2]==1){
    Y1<-Y
  } else if (dim(Y)[2]==2){
    Y1<-Y[,-LeadIndCol]
  } else {
    stop("Y should not have more than 2 columns.")
  }
  if (!is_date_class(est.end.date)){
    stop("est.end.date must be a date class object.")
  }
  if (n.ahead<=0){
    stop("n.ahead must be a positive integer.")
  }
  results <- data.frame(
    Model = names(model_list)
  )
  for (k in 1:n.estimate){
    index_num<-1
    for (model in model_list){
      model$end.date<-est.end.date+(k-1)*gap
      if (class(model)=="SSModelDynamicGompertz"){
        model$Y<-get_timeframe(Y1, model$start.date, model$end.date)
        if (!is.null(model$xpred)){
          model$xpred<-get_timeframe(xpred_targ.full,model$start.date,model$end.date)
        }
        res<-estimate(model)
        if (res$xpred_logical){
          res$xpred.new<-xpred_targ.full
        }
        results[index_num, k+1]=round(mapes(res,n.ahead,Y1)[[criterion]],2)
      } else if (class(model)=="SSModelLeadingIndicator") {
        if (!is.null(model$xpred_lead)){
          model$xpred_lead=xpred_lead.full
        }
        if (!is.null(model$xpred_targ)){
          model$xpred_targ=xpred_targ.full
        }
        res<-estimate(model)
        if (res$xpred_logical[1]){
          res$xpred_lead.new<-xpred_lead.full
        }
        if (res$xpred_logical[2]){
          res$xpred_targ.new<-xpred_targ.full
        }
        results[index_num, k+1]=round(mapes(res,n.ahead,Y)[[criterion]],2)
      } else {
        stop(paste("Model",index_num,"in model_list is not a SSModelDynamicGompertz or SSModelLeadingIndicator object."))
      }
      index_num<-index_num+1
    }
  }
  alldates<-as.character(est.end.date+c(0:(k-1))*gap)
  colnames(results)<-c("Model",alldates)
  return(results)
}

#' @title Identify time resolution of given dates
#'
#' @description This function identifies the time resolution of a vector of 
#' dates, where all dates must come from the same class ("yearmon", "yearqtr", 
#' "date"). 
#'
#' @param dates A vector of dates, length at least 2
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' get_time_resolution(zoo::index(gauteng))
#' 
#' @importFrom timetk is_date_class
#' 
#' @returns A character string from "daily", "monthly", "quarterly" and "yearly".
#'
#' @export
get_time_resolution <- function(dates) {
  if (!is_date_class(dates) || length(dates)<=1){
    stop("Input is not a vector of dates.")
  }
  # Ensure dates are sorted
  dates <- sort(unique(dates))
  date_diff<-diff(dates)  
  if (!isTRUE(all.equal(min(date_diff),max(date_diff)))){
    stop("The dates are not separated by the same time resolution.")
  } else if (class(dates)=="yearqtr"){
    if (min(diff(dates))==1){
      return('yearly')
    } else {
      return('quarterly')
    }
  } else if (class(dates)=="yearmon"){
    if (min(diff(dates))==1){
      return('yearly')
    } else {
      return('monthly')
    }
  } else if (class(dates)=="Date"){
    if (min(diff(dates))==1){
      return('daily')
    }
  } else {
    stop("Input 'dates' must be from classes 'Date', 'yearmon' or 'yearqtr'.")
  }
}

#' @title Change yearqtr object into date object 
#'
#' @description This function changes the yearqtr or yearmon object into a date object, 
#' e.g. the first day of the quarter/month, in order to make plotting these dates 
#' possible in ggplot.
#'
#' @param dates A vector of dates.
#' @returns A vector of transformed dates
#' 
#' @importFrom zoo as.yearmon
#' @examples
#' library(xts)
#' my.date<-zoo::yearqtr(2024.25)
#' qtr2date(my.date)
#' 
#' my.date2<-zoo::yearmon(c(2024.5,2024.8))
#' qtr2date(my.date2)
#'
#' @export
qtr2date<-function(dates){
  return(as.Date(format(as.yearmon(dates), format="%Y-%m-%d")))
}

#' @title Regularly spaced time objects
#'
#' @description Creates a sequence of equally spaced dates that is of the class
#' "yearqtr", "date" or "yearmon". 
#'
#' @param from The first date in the sequence.
#' @param resolution The time resolution of the date sequence. Options are
#' "daily", "monthly", "quarterly" and "yearly".
#' @param to The end date in the sequence, optional. If supplied, \code{to} must be after 
#' (later than) \code{from}. Must be supplied when \code{length.out} is not provided.
#' @param length.out An integer for the length of the sequence, optional. 
#' Must be supplied when \code{to} is not provided. Should not specify both 
#' \code{length.out} and \code{to}.
#' 
#' @importFrom zoo as.yearmon
#' @importFrom zoo as.yearqtr
#' 
#' @returns A vector of dates.
#' @examples
#' #Daily frequency
#' seq_dates(as.Date("2024-01-05"), "daily", length.out=14)
#' 
#' #Quarterly frequency
#' seq_dates(zoo::yearqtr(2020), "quarterly", length.out=12)
#' seq_dates(zoo::yearqtr(2020), "quarterly", to=yearqtr(2022))
#' 
#' #Monthly frequency
#' seq_dates(zoo::yearmon(2020), "monthly", length.out=12)
#'
#' @export
seq_dates<-function(from, resolution, to=NA, length.out=NA){
  if (is.na(length.out) && is.na(to)){
    stop("Both length.out and to inputs cannot be empty.")
  } else if (!is.na(length.out) && !is.na(to)){
    stop("Please supply only one of length.out or to.")
  } else if (is.na(length.out)){
    if (resolution=='daily'){
      seq(from, to, by = 'day')
    } else if (resolution=='quarterly'){
      as.yearqtr(seq(as.numeric(from),
                     as.numeric(to),
                     by=0.25))
    } else if (resolution=='yearly'){
      as.yearmon(seq(as.numeric(from),
                     as.numeric(to),
                     by=1))
    } else if (resolution=='monthly'){
      as.yearmon(seq(as.numeric(from),
                     as.numeric(to),
                     by=1/12))}
  } else if (is.na(to)) {
    if (resolution=='daily'){
      seq(from, by = 'day', length.out = length.out)
    } else if (resolution=='quarterly'){
      as.yearqtr(seq(as.numeric(from), by=0.25, 
                     length.out=length.out))
    } else if (resolution=='yearly'){
      as.yearmon(seq(as.numeric(from), by=1, 
                     length.out=length.out))
    } else if (resolution=='monthly'){
      as.yearmon(seq(as.numeric(from), by=1/12, 
                     length.out=length.out))}
  }
}

