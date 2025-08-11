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
#' @title Class for leading indicator state space model object
#'
#' @description Class for leading indicator state space model object. The model 
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
#' @field xpred1 An xts object containing the values of exogenous variables for 
#' the leading indicator. Dataset must contain values for all dates in the 
#' estimation time frame.
#' @field xpred2 An xts object containing the values of exogenous variables for 
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
    xpred1 = "ANY",  
    xpred2 = "ANY",
    start.date = "ANY",
    end.date = "ANY"),
  methods = list(
    initialize = function(Y, n.lag, sea.period=7, q = NULL,
                          LeadIndCol=1, xpred1=NULL, xpred2=NULL,
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
      xpred1<<-xpred1
      xpred2<<-xpred2
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
      
      y$newCases = lag(y$newCases,n.lag)
      y$LDLcases = lag(y$LDLcases,n.lag)
      y$cCases = lag(y$cCases,n.lag)
      
      y[is.infinite(y)] <- NA
      
      y <- get_timeframe(na.omit(y),start.date)
      
      data_ldl <- get_timeframe(y, start.date, end.date)[,c("LDLcases","LDLhosp")]

      data_mat <- as.matrix(data_ldl)
      
      if (!is.null(xpred1)){
        xpred1<<-get_timeframe(lag(xpred1,n.lag),index(data_ldl)[1],tail(index(data_ldl),1))
      }
      if (!is.null(xpred2)){
        xpred2<<-get_timeframe(xpred2,index(data_ldl)[1],tail(index(data_ldl),1))
      }
      
      # update = function(pars, model, q) {
      #   "Update method for Kalman filter to implement the dynamic Gompertz curve
      #  model.
      #  A maximum of 3 parameters are used to set the observation noise
      #  (1 parameter), the transition equation slope and seasonal noise. If q (signal
      #   to noise ratio) is not null then the slope noise is set using this
      #   ratio.
      #  \\subsection{Parameters}{\\itemize{
      #   \\item{\\code{pars} Vector of parameters.}
      #   \\item{\\code{model} \\code{KFS} model object.}
      #   \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
      #    variance).}
      # }}
      # \\subsection{Return Value}{\\code{KFS} model object.}"
      #   estH <- any(is.na(model$H))
      #   estQ <- any(is.na(model$Q))
      #   if ((!estH) & (!estQ)) {
      #     # If nothing to update then return model
      #     return(model)
      #   } else {
      #     nparQ <- 0
      #     # 1. Set seasonal noise
      #     if (estQ) {
      #       Q <- as.matrix(model$Q[, , 1])
      #       T <- as.matrix(model$T[, , 1])
      #       
      #       # Update diagonal elements
      #       naQd <- which(is.na(diag(Q)))
      #       Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
      #       
      #       # Check for off-diagonal elements and raise error if found.
      #       naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
      #       if (length(naQnd) > 0) {
      #         stop("NotImplmentedError: Unexpected off-diagonal element updating")
      #       }
      #       
      #       if (sea.period >1){
      #         #Identify elements corresponding to seasonal cases 
      #         nparQ <- 1
      #         na_sea_cases<-grep("^sea\\_trig.*LDLcases$",rownames(T))
      #         diag(Q)[na_sea_cases] <- exp(0.5 * pars[nparQ])
      #         #Identify elements corresponding to seasonal hospitalizations
      #         nparQ <- nparQ+1
      #         na_sea_hosp<-grep("^sea\\_trig.*LDLhosp$",rownames(T))
      #         diag(Q)[na_sea_hosp] <- exp(0.5 * pars[nparQ])
      #       }
      #       
      #       # 2. Set observation noise
      #       H <- as.matrix(model$H[, , 1])
      #       if (estH) {
      #         naHd <- which(is.na(diag(H)))
      #         H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
      #         nparQ<-nparQ+1
      #         diag(H)[naHd] <- exp(0.5 * pars[nparQ])
      #         model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
      #       }
      #       
      #       # 3. Set slope noise
      #       # Get index of slope
      #       model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
      #       i.slope <- grep("slope",rownames(T))
      #       # Estimate slope if no signal to noise ratio specified.
      #       if (is.null(q)) {
      #         nparQ<-nparQ+1
      #         Q.slope <- exp(0.5 * pars[nparQ])
      #       } else {
      #         Q.slope <- crossprod(H[naHd, naHd]) * q
      #       }
      #       model$Q[i.slope, i.slope, 1] <- Q.slope
      #       
      #       # 4. Set AR1 noise
      #       if (ar1){
      #         nparQ<-nparQ+1
      #         i.ar1 <- nrow(Q)
      #         Q[i.ar1, i.ar1] <- exp(0.5 * pars[nparQ])
      #         model$Q[i.ar1, i.ar1, 1] <- Q[i.ar1, i.ar1]
      #         
      #         nparQ<-nparQ+1
      #         T <- model$T[,,1]
      #         model$T[nrow(T),ncol(T),1] <- pars[nparQ]
      #       }
      #     }
      #   }
      #   return(model)
      # }

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
      # an extra trend [random walk] in LDLhosp only [degree = 1],
      # and 7 day dummy variable seasonal.
      
      if (sea.period<2){
        if (is.null(xpred1)){
          if (is.null(xpred2)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred2, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        } else {
          if (is.null(xpred2)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred1, type="distinct", index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred1, type="distinct", index=1)+
                             SSMregression(~xpred2, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        }
      }
      else {
        if (is.null(xpred1)){
          if (is.null(xpred2)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred2, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        } else {
          if (is.null(xpred2)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred1, type="distinct", index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xpred1, type="distinct", index=1)+
                             SSMregression(~xpred2, type="distinct", index=2),
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
        update = updatesn %>% partial(snr=q,order=2,index=1)

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
        xpred_logical=c(!is.null(xpred1),!is.null(xpred2)),
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
      #   qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
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
        geom_line(aes(y = newCases, color = series.name.lead), lwd = 0.85) +
        geom_line(aes(y = newAdmit, color = series.name.target), lwd = 0.85) +
        scale_color_manual(values = c("red", "blue"))+
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.2, 0.85),
          legend.title = element_text(size = 5),
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
