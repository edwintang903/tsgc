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
#' @title  Class for leading indicator state space model object
#'
#' @description Class for leading indicator state space model object.
#'
#' \subsection{Methods}{
#' \code{get_model(y, q = NULL, sea.type = 'trigonometric', sea.period = 7)}
#' Retrieves the model object.
#' \subsection{Parameters}{\itemize{
#'  \item{\code{y} The cumulated variable.}
#'  \item{\code{q} The signal-to-noise ratio (ratio of slope to irregular
#' variance). Defaults to \code{'NULL'}, in which case no signal-to-noise ratio
#' will be imposed. Instead, it will be estimated.}
#'  \item{\code{sea.type} Seasonal type. Options are \code{'trigonometric'} and
#' \code{'none'}. \code{'trigonometric'} will yield a model with a trigonometric
#' seasonal component and \code{'none'} will yield a model with no seasonal
#' component.}
#'  \item{\code{sea.period}  The period of seasonality. For a day-of-the-week
#' effect with daily data, this would be 7. Not required if
#' \code{sea.type = 'none'}.}
#' }}
#' \subsection{Return Value}{\code{KFS} model object.}
#' }
#'
#' @importFrom xts periodicity last
#' @importFrom methods new
#' @importFrom magrittr %>%
#' @importFrom KFAS SSMtrend SSMseasonal SSModel
#' @importFrom purrr partial
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
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
    LeadIndCol ="numeric"
  ),
  methods = list(
    initialize = function(Y, n.lag,sea.period=7, q = NA,LeadIndCol=1)
    {
      "Create an instance of the \\code{SSModelLeadingIndicator} class.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{Y} The cumulated leading indiactor variable and target variable.}
        \\item{\\code{n.lag} Number of days to lag the leading indicator.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance). Defaults to \\code{'NULL'}, in which case no
         signal-to-noise ratio will be imposed. Instead, it will be estimated.}
         \\item{\\code{LeadIndCol} The column in \\code{Y} that contains the leading indicator}
      }}
      \\subsection{Usage}{\\code{SSModelLeadingIndicator(Y=y,n.lag = n.lag,q=q,LeadIndCol=1)}}"
      Y <<- Y
      q <<- q
      sea.period<<-sea.period
      n.lag <<- n.lag
      LeadIndCol <<- LeadIndCol
    },
    estimate = function()
    {
    "Retrieves the model object."

      # Compute LDL and lag data appropriately
      y<-add_daily_ldl(Y,LeadIndCol=LeadIndCol)

      data_ldl = y[,c("LDLcases","LDLhosp")]

      data_ldl$LDLcases = lag(as.vector(data_ldl$LDLcases),n.lag)

      data_ldl <- na.omit(data_ldl)

      data_mat = as.matrix(data_ldl)

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

      if (is.na(sea.period)){
        mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                      SSMtrend(degree = 1, Q = matrix(NA),index=1),
                      H = matrix(c(NA,0,0,NA),2,2))
      }
      else {
        mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                      SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='dummy', type='distinct')+
                      SSMtrend(degree = 1, Q = matrix(NA),index=1),
                      H = matrix(c(NA,0,0,NA),2,2))
      }

      # Compute number of parameters - this is just the number of NAs in the model Q and H combined.
      npar = sum(is.na(mod$Q)) + sum(is.na(mod$H))

      # Set the options for the update function
      # We have a signal/noise ratio of 0.005, the signal is the slope and we are
      # targeting the variance of the irregular in cases

      if (is.na(q)){
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
        index= index(y),
        data_xts = y,
        output = out,
        n.lag=n.lag,
        sea.period=sea.period,
        LeadIndCol=LeadIndCol)
      return(results)},
    summary = function() {
      out <- output(.self$estimate())
      # q<-.self$q
      # if(is.null(q)){
      #   qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      # }
      dates<-index(.self$Y)
      
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
      cat("  - Dataset start date:", format(as.Date(dates[1], origin = "1970-01-01")))
      cat("\n")
      cat("  - Dataset end date:", format(as.Date(tail(dates,1), origin = "1970-01-01")))
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(out)
    },
    print = function() {
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
    plot_diff =function(...){
      plot(diff(.self$Y), ...)
    }
  )
)
