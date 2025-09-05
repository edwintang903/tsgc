# Created by: Craig Thamotheram
# Created on: 27/07/2022

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

setOldClass("KFS")
#'
#' @title Class for designing a Leading Indicator Model
#'
#' @description A class for specifying the parameters of a leading indicator model. The model 
#' settings are stored in the fields of this object, and the class contains 
#' methods to obtain FilterResultsLI object for further analysis and plotting 
#' the time series under investigation as an exploratory data analysis.
#'
#' @field Y The cumulated variable.
#' @field q The signal-to-noise ratio (ratio of slope to irregular
#'   variance). Defaults to \code{'NULL'}, in which case no
#'   signal-to-noise ratio will be imposed. Instead, it will be estimated.
#'@field sea.period The period of seasonality. For a day-of-the-week
#'   effect with daily data, this would be 7. Not required if
#'   \code{sea.type = 'none'}.
#' @field n.lag Number of days/months/quarters/years to lag the leading indicator.
#' @field xpred_lead An xts object containing the values of exogenous variables for 
#' the leading indicator. Dataset must contain values for all dates in the 
#' estimation time frame.
#' @field xpred_targ An xts object containing the values of exogenous variables for 
#' the target variable. Dataset must contain values for all dates in the 
#' estimation time frame.
#' @field LeadIndCol The column in \code{Y} that contains the leading indicator.
#' @field start.date Start date of the estimation period for estimating the target variable. 
#' Must be one of the following types: \code{yearqtr}, \code{date} or \code{yearmon}. 
#' @field end.date End date of the estimation period for estimating the target variable. 
#' Must be one of the following types: \code{yearqtr}, \code{date} or \code{yearmon}. 
#'
#' @importFrom xts periodicity last lag.xts
#' @importFrom methods new
#' @importFrom magrittr %>%
#' @importFrom KFAS SSMtrend SSMseasonal SSModel
#' @importFrom purrr partial
#'
#' @examples
#' library(tsgc)
#' 
#' # Specify a model with the estimation timeframe
#' estimation.date.start <- as.Date("2020-02-25")
#' estimation.date.end <- as.Date("2020-04-01")
#' 
#' model <- SSModelLeadingIndicator(Y=ukitaly, n.lag = 14, 
#' sea.period = 7,LeadIndCol=1, start.date=estimation.date.start,
#' end.date=estimation.date.end)
#' 
#' # Show summary of the model object
#' summary(model)
#' 
#' # Print a short description of the model object
#' print(model)
#' 
#' # Plot the time series stored in the model object
#' plot(model, title="COVID Daily Cases in UK and Italy",
#' series.name.lead="Italy", series.name.target="UK", take.log=FALSE)
#'
#' # Estimate a specified model
#' res <- estimate(model)
#' res
#'
#' @export SSModelLeadingIndicator
#' @exportClass SSModelLeadingIndicator
SSModelLeadingIndicator <- setRefClass(
  "SSModelLeadingIndicator",
  fields = list(
    Y = "ANY",
    q = "ANY",
    sea.period= "numeric",
    n.lag = "numeric",
    LeadIndCol ="numeric",
    xpred_lead = "ANY",  
    xpred_targ = "ANY",
    start.date = "ANY",
    end.date = "ANY"),
  methods = list(
    initialize = function(Y, n.lag, sea.period=7, q = NULL,
                          LeadIndCol=1, xpred_lead=NULL, xpred_targ=NULL,
                          start.date=index(Y)[1], end.date=tail(index(Y),1))
    {
      "Create an instance of the \\code{SSModelLeadingIndicator} class with the 
      fields laid out at the beginning of the documentation."
      resolu<-get_time_resolution(index(Y))
      Y <<- Y
      q <<- q
      sea.period<<-sea.period
      n.lag <<- if (resolu=="daily" || resolu=="yearly"){
        n.lag
      } else if (resolu=="quarterly"){
        n.lag*4
      } else if (resolu=="monthly"){
        n.lag*12}
      LeadIndCol <<- LeadIndCol
      xpred_lead<<-xpred_lead
      xpred_targ<<-xpred_targ
      start.date<<-start.date
      end.date<<-end.date
    },
    estimate = function()
    {
    "Estimates the Leading Indicator model when applied to an object of
      class \\code{SSModelLeadingIndicator}.
      \\subsection{Return Value}{An object of class \\code{FilterResultsLI}
      containing the result output for the estimated Leading Indicator
      model.}"
      
      # Compute LDL and lag data appropriately
      y<-add_daily_ldl(Y, LeadIndCol=LeadIndCol)
      
      y$newLead = stats::lag(y$newLead,n.lag)
      y$LDLlead = stats::lag(y$LDLlead,n.lag)
      y$cLead = stats::lag(y$cLead,n.lag)
      
      y[is.infinite(y)] <- NA
      
      y <- get_timeframe(na.omit(y),start.date)
      
      data_ldl <- get_timeframe(y, start.date, end.date)[,c("LDLlead","LDLtarg")]

      data_mat <- as.matrix(data_ldl)
      
      if (!is.null(xpred_lead)){
        xpred_lead<<-get_timeframe(stats::lag(xpred_lead,n.lag),index(data_ldl)[1],tail(index(data_ldl),1))
      }
      if (!is.null(xpred_targ)){
        xpred_targ<<-get_timeframe(xpred_targ,index(data_ldl)[1],tail(index(data_ldl),1))
      }

      # Standard update function - edited to allow the targeting of the signal-to-noise ratio
      # Signal-to-noise ratio is defined as the variance of the trend component of order 'order'
      # (= 1 for level, = 2 for slope, etc) relative to variance of irregular of series 'index'
      # (= 1 for 1st col of dataframe, = 2 for 2nd etc)
      updatesn=function(pars, model, snr, order, index){
        if(any(is.na(model$Q))){
          Q <- as.matrix(model$Q[,,1])
          naQd  <- which(is.na(diag(Q)))
          naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))
          Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
          
          diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
          Q[naQd,naQd][naQnd] <- pars[length(naQd)+1:length(naQnd)]
          model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
        }
        if(!identical(model$H,'Omitted') && any(is.na(model$H))){
          H<-as.matrix(model$H[,,1])
          naHd  <- which(is.na(diag(H)))
          naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
          H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0
          diag(H)[naHd] <-
            exp(0.5 * pars[length(naQd)+length(naQnd)+1:length(naHd)])
          H[naHd,naHd][naHnd] <-
            pars[length(naQd)+length(naQnd)+length(naHd)+1:length(naHnd)]
          model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
          model$Q[order,order,1] <- snr*crossprod(H[index,index])
        }
        model
      }
      # Create the SSM model
      # This has a common trend and slope (common trend of degree 2),
      # an extra trend [random walk] in LDLtarg only [degree = 1],
      # and 7 day dummy variable seasonal.
      
      if (sea.period<2){
        if (is.null(xpred_lead)){
          if (is.null(xpred_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        } else {
          if (is.null(xpred_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred_lead, type="distinct", index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred_lead, type="distinct", index=1)+
                             SSMregression(~xpred_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        }
      }
      else {
        if (is.null(xpred_lead)){
          if (is.null(xpred_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        } else {
          if (is.null(xpred_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred_lead, type="distinct", index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred_lead, type="distinct", index=1)+
                             SSMregression(~xpred_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        }
      }

      # Compute number of parameters - this is just the number of NAs in the model Q and H combined.
      npar = sum(is.na(mod$Q)) + sum(is.na(mod$H))

      # Set the options for the update function
      # We have a signal/noise ratio of 0.005, the signal is the slope and we are
      # targeting the variance of the irregular in cases

      if (is.null(q)){
        fit = fitSSM(mod, rep(0,npar))
      }
      else{
        update = updatesn %>% partial(snr=q,order=2,index=2)

        # Fit the state-space model (ML, diffuse prior)
        fit = fitSSM(mod, rep(0,npar), updatefn = update)
      }

      # Apply the Kalman filter and smoother to the fitted model
      out = KFS(fit$model)

      results <- FilterResultsLI$new(
        data_xts = y,
        output = out,
        n.lag=n.lag,
        sea.period=sea.period,
        LeadIndCol=LeadIndCol,
        xpred_logical=c(!is.null(xpred_lead),!is.null(xpred_targ)),
        start.date=index(data_ldl)[1],
        end.date=tail(index(data_ldl),1))
      return(results)},
    summary = function() {
      "Supplies details of the SSModelLeadingIndicator object, such as estimated 
      parameter values, start and end dates of estimation."
      result<-.self$estimate()
      out <- output(result)
      # q<-.self$q
      # if(is.null(q)){
      #   qest <- matrixKFS(out,"Q")[3, 3, 1]/matrixKFS(out,"H")[2, 2, 1]
      # }
      start<-result$start.date
      end<-result$end.date
      resolution<-result$resolution
      
      cat("Summary of SSModelLeadingIndicator Model")
      cat("\n")
      cat("--------------------------------------\n")
      cat("Cumulated Variable:\n")
      base::print(head(.self$Y))
      # cat("Signal-to-Noise Ratio (q):",
      #     ifelse(is.null(q), paste(signif(qest,3), "(estimated)"),
      #            paste(q, ("(user specified)"))), "\n")
      cat("Model Details:\n")
      cat("  - Model Type: Leading Indicator Model")
      cat("\n")
      cat("  - Seasonal Component: ", ifelse(is.na(sea.period), "None", "Trigonometric"), "\n")
      cat("  - Period of Seasonality: ", ifelse(is.na(sea.period), "N/A", sea.period), "\n")
      if (resolution=="daily"){
        cat("  - Estimation start date:", format(as.Date(start, origin = "1970-01-01"))) 
        cat("\n")
        cat("  - Estimation end date:", format(as.Date(end, origin = "1970-01-01")))
      } else if (resolution=="quarterly"){
        cat("  - Estimation start date:", format(as.yearqtr(start))) 
        cat("\n")
        cat("  - Estimation end date:", format(as.yearqtr(end)))
      } else if (resolution=="monthly"  || resolution=="yearly"){
        cat("  - Estimation start date:", format(as.yearmon(start))) 
        cat("\n")
        cat("  - Estimation end date:", format(as.yearmon(end)))
      } 
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(out)
    },
    print = function() {
      "Provides a quick description of SSModelDynamicGompertz object, providing 
      model states and standard errors."
      
      out <- output(.self$estimate()) #KFS object
      # if(is.null(.self$q)){
      #   qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      # }
      cat("SSModelLeadingIndicator Model")
      cat("\n")
      cat("\n")
      cat("Cumulated Variable:\n")
      base::print(head(.self$Y))
      cat("Number of observations:", length(.self$Y))
      cat("\n")
      # cat("Signal-to-Noise Ratio (q):", 
      #     ifelse(is.null(.self$q), paste(signif(qest,5), "(estimated)"), 
      #            paste(.self$q, ("(user specified)"))), "\n")
      cat("Seasonal components?",
          ifelse(is.null(seasonalComp(out)),
                 "No","Yes"),"\n")
    },
    plot=function(title=NULL, series.name.lead="Leading Indicator", 
                  series.name.target="Target Variable",
                  date_break=NULL, take.log=TRUE){
      "Plots the lagged differences of the cumulated dataset \\code{Y} in this 
      \\code{SSModelLinearIndicator} object against time, which could represent 
      daily cases.
      \\subsection{Parameters}{\\itemize{
        \\item{\\code{title} Title for forecast plot. Enter as character string. 
        \\code{NULL} (i.e. no title) by default.}
        \\item{\\code{series.name.lead} The name of the leading indicator series 
        for. E.g. \\code{'cases'}. Enter as character string. Default is `Leading Indicator`.}
        \\item{\\code{series.name.target} The name of the target variable series 
        for. E.g. \\code{'hospitalizations'}. Enter as character string. Default is `Target Variable`.}
        \\item{\\code{date_break} A character string (e.g. '60 days') specifying the interval 
        between date labels, used in \\code{scale_x_date} within 
        \\code{ggplot}. If \\code{NULL} (default), a suitable interval is chosen 
        automatically by \\code{ggplot}.}
        \\item{\\code{take.log} A logical value indicating whether to return 
        take log of the lagged differences. Defaults to \\code{TRUE}.}
        }
      }
      \\subsection{Return Value}{A plot of the lagged differences of the 
      cumulated dataset \\code{Y} against time.}"
      # Transform the data to calculate daily cases and log growth rates.
      eng_full <- add_daily_ldl(Y)
      eng_daily <- eng_full[, 3:4]
      dates<-index(eng_daily)
      reso<-get_time_resolution(dates)
      if (reso!="daily"){
        dates<-qtr2date(dates)
      }
      
      # Plot daily new cases and admissions.
      if (take.log){
        base_plot<-ggplot(log(eng_daily), aes(x = dates))+
          labs(
            title = title,
            x = "Date",
            y = "log(Number)",
            color = "Legend"
          ) 
      } else {
        base_plot<-ggplot(eng_daily, aes(x = dates))+
          labs(
            title = title,
            x = "Date",
            y = "Number",
            color = "Legend"
          ) 
      }
      
      data_plot<-base_plot+
        geom_line(aes(y = newLead, color = series.name.lead), lwd = 0.85) +
        geom_line(aes(y = newTarg, color = series.name.target), lwd = 0.85) +
        scale_color_manual(values = c("red", "blue"))+
        theme(
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold")
        ) 
      if (!is.null(date_break)) {
        data_plot <- data_plot + scale_x_date(date_breaks = date_break)
      } 
      data_plot
    }
  )
)
