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

setOldClass("xts")
setOldClass("KFS")
#'
#' @title  Class for dynamic Gompertz curve state space model object.
#'
#' @description Class for dynamic Gompertz curve state space model object.
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
#' @importFrom xts xts
#' @importFrom zoo index
#' @importFrom KFAS SSModel fitSSM KFS
#' @importFrom magrittr %>%
#' @importFrom KFAS SSMtrend SSMseasonal
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
#' @export SSModelDynamicGompertz
#' @exportClass SSModelDynamicGompertz
SSModelDynamicGompertz <- setRefClass(
  "SSModelDynamicGompertz",
  fields = list(
    Y = "xts",
    q = "ANY",  # No native option for numeric | NULL - see
    # https://stackoverflow.com/questions/24363069/multiple-acceptable-classes-
    # in-reference-class-field-list
    sea.type="ANY", 
    sea.period="ANY",
    reinit.date = "ANY",
    original.results = "ANY",
    use.presample.info = "ANY"
  ),
  methods = list(initialize = function(Y, q = NULL, sea.type = 'trigonometric',
                                       sea.period = 7,reinit.date=NULL, original.results=NULL,
                                       use.presample.info=TRUE)
  {
    "Create an instance of the \\code{SSModelDynGompertzReinit} class.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{Y} The cumulated variable.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance). Defaults to \\code{'NULL'}, in which case no
         signal-to-noise ratio will be imposed. Instead, it will be estimated.}
        \\item{\\code{reinit.date} The reinitialisation date \\eqn{r}. Should
        be specified as an object of class \\code{\"Date\"}. Must be specified
        if \\code{original.results = NULL} and
        \\code{use.pre.sample.info = TRUE}.}
        \\item{\\code{original.results} In place of a reinitialisation date, a
        \\code{KFS} results object can be specified here and the parameters for
         the reinitialisation will be taken from this object. Must be specified
          if \\code{reinit.date = NULL} and \\code{use.pre.sample.info = TRUE}.}
        \\item{\\code{use.presample.info}  Logical value denoting whether or not
         to use information from before the reinitialisation date in the
         reinitialisation procedure. Default is \\code{TRUE}. If \\code{FALSE},
         the model is estimated from scratch from the reinitialisation date and
         no attempt to use information from before the reinitialisation date is
          made.}
      }}
      \\subsection{Usage}{\\code{SSModelDynGompertzReinit$new(y, q = 0.005,
      reinit.date = as.Date(\"2021-05-12\",format = date.format))}}"
    Y <<- Y
    q <<- q
    sea.type <<- sea.type
    sea.period <<- sea.period
    reinit.date <<- reinit.date
    original.results <<- original.results
    use.presample.info <<- use.presample.info
  },get_dynamic_gompertz_model = function(
    y,
    q = NULL,
    sea.type = 'trigonometric',
    sea.period = 7,
    a1 = NULL,
    P1 = NULL,
    Q = NULL,
    H = NULL
  )
  {
    "Returns dynamic Gompertz curve model.
    \\subsection{Parameters}{\\itemize{
      \\item{\\code{y} The cumulated variable}
      \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
      variance). Defaults to \\code{'NULL'}, in which case no signal-to-noise
      ratio will be imposed. Instead, it will be estimated.}
      \\item{\\code{sea.type} Seasonal type. Options are \\code{'trigonometric'}
       and \\code{'none'}. \\code{'trigonometric'} will yield a model with a
       trigonometric seasonal component and \\code{'none'} will yield a model
       with no seasonal component.}
      \\item{\\code{sea.period} The period of seasonality. For a day-of-the-week
       effect with daily data, this would be 7. Not required if
       \\code{sea.type = 'none'}.}
      \\item{\\code{a1} Optional parameter specifying the prior mean of the
      states. Defaults to \\code{'NULL'}. Leave as \\code{'NULL'} for a diffuse
      prior (no prior information). If a proper prior is to be specified, both
      \\code{a1} and \\code{P1} must be given.}
      \\item{\\code{P1} Optional parameter specifying the prior mean of the
      states. Defaults to \\code{'NULL'}. Leave as \\code{'NULL'} for a diffuse
       prior (no prior information). If a proper prior is to be specified,
       both \\code{a1} and \\code{P1} must be given.}
      \\item{\\code{Q} Optional parameter specifying the state error variances
      where these are to be imposed rather than estimated. Defaults to
      \\code{'NULL'} which will see the variances estimated.}
      \\item{\\code{H} Optional parameter specifying the irregular variance
      where this is to be imposed rather than estimated. Defaults to
      \\code{'NULL'} which will see the variance estimated.}
    }}
    \\subsection{Description}{
    The dynamic Gompertz with an integrated random walk (IRW) trend is
    \\deqn{\\ln g_{t}=\\delta_{t}+\\varepsilon_{t},  \\;\\;\\;\\;
    \\varepsilon_{t}\\sim NID(0,\\sigma_{\\varepsilon }^{2}), \\;\\;\\;\\;
    t=2,...,T, }
    where \\eqn{Y_t} is the cumulated variable, \\eqn{y_t = \\Delta Y_t},
    \\eqn{\\ln g_{t}=\\ln y_{t}-\\ln Y_{t-1}} and
    \\deqn{\\delta_{t} =\\delta_{t-1}+\\gamma_{t-1},}
    \\deqn{\\gamma_{t} =\\gamma_{t-1}+\\zeta_{t}, \\;\\;\\;\\;
    \\zeta_{t}\\sim NID(0,\\sigma_{\\zeta }^{2}),}
    where the observation disturbances \\eqn{\\varepsilon_{t}}  and slope
    disturbances \\eqn{\\zeta_{t}}, are iid Normal and mutually independent.
    Note that, the larger the signal-to-noise ratio,
    \\eqn{q_{\\zeta }=\\sigma_{\\zeta }^{2}/\\sigma_{\\varepsilon }^{2}},
    the faster the slope changes in response to new observations. Conversely,
    a lower signal-to-noise ratio induces smoothness.

    For the model without seasonal terms (\\code{sea.type = 'none'}) the are
    priors are
    \\deqn{\\begin{pmatrix} \\delta_1 \\ \\gamma_1 \\end{pmatrix}
    \\sim N(a_1,P_1)}.
    The diffuse prior has \\eqn{P_1 = \\kappa I_{2\\times 2}} with
    \\eqn{\\kappa \\to \\infty}. Implementation of the diffuse prior is handled
     by the package \\code{KFAS} (Helske, 2017). Where the model has a seasonal
      component (\\code{sea.type = 'trigonometric'}), the vector of prior means
       \\eqn{a_1} and the prior covariance matrix \\eqn{P_1} are extended
       accordingly.

    See the vignette for details of the variance matrix \\eqn{Q}.
    \\eqn{H = \\sigma^2_{\\varepsilon}}.
    }
    "
    Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
    Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
    Ht <- if (is.null(H)) { NA } else { H }
    
    # 1. Set prior on state as ~ N(a1, P1) if a1 supplied.
    use.prior <- if (!is.null(a1)) { TRUE } else { FALSE }
    
    if (use.prior) {
      if (sea.type == 'trigonometric') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope)),
              a1 = a1[1:2],
              P1 = P1[1:2, 1:2]
            ) +
            SSMseasonal(
              period = sea.period,
              Q = Qt.seas,
              sea.type = sea.type,
              a1 = a1[3:dim(a1)[1]],
              P1 = P1[3:dim(a1)[1], 3:dim(a1)[1]]
            ),
          H = Ht
        )
        n.pars <- 0
      } else if (sea.type == 'none') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope)),
              a1 = a1[1:2],
              P1 = P1[1:2, 1:2]
            ),
          H = Ht
        )
        n.pars <- 0
      } else {
        stop(sprintf("sea.type= '%s' not implemented", sea.type))
      }
    } else {
      if (sea.type == 'trigonometric') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope))
            ) +
            SSMseasonal(
              period = sea.period,
              Q = Qt.seas,
              sea.type = sea.type),
          H = matrix(Ht)
        )
        n.pars <- if (is.null(q)) { 3 } else { 2 }
      } else if (sea.type == 'none') {
        ss_model <- SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = list(matrix(0), matrix(Qt.slope))
            ),
          H = matrix(Ht)
        )
        n.pars <- if (is.null(q)) { 2 } else { 1 }
      } else {
        stop(sprintf("sea.type= '%s' not implemented", sea.type))
      }
    }
    out <- list(model = ss_model, inits = rep(0, n.pars))
    return(out)
  },
  update = function(pars, model, q, sea.type) {
    "Update method for Kalman filter to implement the dynamic Gompertz curve
       model.
       A maximum of 3 parameters are used to set the observation noise
       (1 parameter), the transition equation slope and seasonal noise. If q (signal
        to noise ratio) is not null then the slope noise is set using this
        ratio.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{pars} Vector of parameters.}
        \\item{\\code{model} \\code{KFS} model object.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance).}
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}.}
      }}
      \\subsection{Return Value}{\\code{KFS} model object.}"
    estH <- any(is.na(model$H))
    estQ <- any(is.na(model$Q))
    if ((!estH) & (!estQ)) {
      # If nothing to update then return model
      return(model)
    } else {
      nparQ <- if (sea.type == 'trigonometric') { 1 } else { 0 }
      # 1. Set seasonal noise
      if (estQ) {
        Q <- as.matrix(model$Q[, , 1])
        # Update diagonal elements
        naQd <- which(is.na(diag(Q)))
        Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
        diag(Q)[naQd] <- exp(0.5 * pars[1])
        # Check for off-diagonal elements and raise error if found.
        naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
        if (length(naQnd) > 0) {
          stop("NotImplmentedError: Unexpected off-diaganol element updating")
        }
      }
      
      # 2. Set observation noise
      H <- as.matrix(model$H[, , 1])
      if (estH) {
        naHd <- which(is.na(diag(H)))
        H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
        diag(H)[naHd] <- exp(0.5 * pars[(nparQ + 1)])
        model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
      }
      
      # 3. Set slope noise
      # Get index of slope, 1 before the seasonal component.
      model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
      i.slope <- 2
      # Estimate slope if no signal to noise ratio specified.
      if (is.null(q)) {
        Q.slope <- exp(0.5 * pars[(nparQ + 2)])
      } else {
        Q.slope <- crossprod(H[naHd, naHd]) * q
      }
      model$Q[i.slope, i.slope, 1] <- Q.slope
    }
    return(model)
  },
  estimate = function() {
    "Estimates the dynamic Gompertz curve model when applied to an object of
      class \\code{SSModelDynamicGompertz} or \\code{SSModelDynGompertzReinit}.
      \\subsection{Parameters}{\\itemize{
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}. \\code{'trigonometric'} will
         yield a model with a trigonometric seasonal component and
         \\code{'none'} will yield a model with no seasonal component.}
        \\item{\\code{sea.period} The period of seasonality. For a
        day-of-the-week effect with daily data, this would be 7. Not required
        if \\code{sea.type = 'none'}.}
      }}
      \\subsection{Return Value}{An object of class \\code{FilterResults}
      containing the result output for the estimated dynamic Gompertz curve
      model.}
      "
    # 1. Get LDL of cumulative series Y.
    y <- tsgc::df2ldl(Y)
    
    # 2. Add update / model methods
    updatefn <- purrr::partial(
      .self$update, ... =, q = q, sea.type = sea.type
    )
    model <- .self$get_model(y, q = q, sea.type, sea.period)
    # 2. Estimate via MLE unknown params
    model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                        method = 'BFGS')
    
    # 3. Run smoother/filter
    model_output <- KFS(model_fit$model)
    
    # 4. Get truncated index from model if using a reinitialisation in
    # self$get_model
    date.index <- if (!is.null(model$index)) { model$index } else { index(y) }
    
    results <- FilterResults$new(
      data_xts = Y,
      index = date.index,
      reinit.date=reinit.date,
      output = model_output
    )
    return(results)
  },
    get_model = function(
    y,
    q = NULL,
    sea.type = 'trigonometric',
    sea.period = 7
    )
    {
      "Retrieves the model object.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{y} The cumulated variable.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
        variance). Defaults to \\code{'NULL'}, in which case no signal-to-noise
        ratio will be imposed. Instead, it will be estimated.}
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}. \\code{'trigonometric'} will
         yield a model with a trigonometric seasonal component and
         \\code{'none'} will yield a model with no seasonal component.}
        \\item{\\code{sea.period}  The period of seasonality. For a
        day-of-the-week effect with daily data, this would be 7. Not required
        if \\code{sea.type = 'none'}.}
      }}
      \\subsection{Return Value}{\\code{KFS} model object.}"
      if (is.null(.self$reinit.date)){
        model <- .self$get_dynamic_gompertz_model(
          y, q = q, sea.type = sea.type, sea.period = sea.period
        )
        return(model)
      }
      else{
        # 4.1. Index for reinitialisation, t_0
        stopifnot(length(Y[reinit.date]) == 1)
        Y.t.r_0 <- as.numeric(Y[reinit.date - 1])
        
        # 4.2 Reinitialisation:
        #   ln g_t^r = ln g_t + ln (Y_{t-1}/Y_{t-1}^r), where Y_t^r=Y_t-Y_{r_0}.
        idx.dates <- (index(y) > reinit.date)
        lag.Y <- stats::lag(Y)[idx.dates]
        y.reinit <- y[index(y) > reinit.date] + log(lag.Y / (lag.Y - Y.t.r_0))
        
        # 4.3 Run Kalman filter/smoother on new series with non-diffuse prior
        if (use.presample.info) {
          # Either estimate full model here or take results from previous model.
          if (is.null(original.results)) {
            # NB. Restrict sample to t<=r - date of reinitialisation.
            idx.est <- zoo::index(Y) <= reinit.date
            model <- SSModelDynamicGompertz$new(Y = Y[idx.est], q = q)
            res.original <- model$estimate()
            model_output <- output(res.original)
          } else {
            model_output <- output(original.results)
            model_seasonal <- seasonalComp(original.results)
            sea.type <<- if (is.null(model_seasonal)) {'none'} else {
              'trigonometric'}
            sea.period <<- if (!is.null(model_seasonal)) {
              ncol(att(model_output)) - 1}
          }
          
          # 4.3 Reset slope to 0 and add constant to initial value for level.
          # where reinit.date is t=r
          idx <- which(reinit.date == index(y))
          stopifnot(length(idx) == 1)
          att <- att(model_output)[idx,]
          Ptt <- Ptt(model_output)[, , idx]
          Tt <- drop(matrixKFS(model_output,"T"))
          Rt <- drop(matrixKFS(model_output,"R"))
          Qt <- drop(matrixKFS(model_output,"Q"))
          Ht <- drop(matrixKFS(model_output,"H"))
          
          # a. Take a_{r|r} and P_{r|r} through prediction step to get a_{r+1}
          # and P_{r+1}
          a1 <- Tt %*% att
          P1 <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
          
          # b. Set slope to 0 and add correction (\ln(Y_r/y_r) to level.
          a1["slope",] <- 0
          a1["level",] <- a1["level",] + log(Y[idx] / (Y[idx] - Y.t.r_0))
        } else {
          # Don't use presample info
          a1 <- NULL; P1 <- NULL; Qt <- NULL; Ht <- NULL
        }
        out <- .self$get_dynamic_gompertz_model(
          y = y.reinit, q = q, sea.type = sea.type, sea.period = sea.period,
          a1 = a1, P1 = P1, Q = Qt, H = Ht
        )
        out[['index']] <- index(y.reinit)
        return(out)
      }
    },
    summary = function() {
      out <- suppressWarnings(output(.self$estimate()))
      q <<-.self$q
      if(is.null(q)){
        qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      }
      reinit<-!is.null(reinit.date)
      dates<-index(Y)
      
      cat("Summary of SSModelDynamicGompertz Model")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("--------------------------------------\n")
      cat("Cumulated Variable:\n")
      base::print(head(.self$Y))
      cat("Signal-to-Noise Ratio (q):", 
          ifelse(is.null(q), paste(signif(qest,3), "(estimated)"), 
                 paste(q, ("(user specified)"))), "\n")
      cat("Model Details:\n")
      cat("  - Model Type: Dynamic Gompertz Curve")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("  - Seasonal Component: ", ifelse(sea.type == 'none', "None", "Trigonometric"), "\n")
      cat("  - Period of Seasonality: ", ifelse(sea.type == 'none', "N/A", sea.period), "\n")
      cat("  - Dataset start date:", format(as.Date(dates[1], origin = "1970-01-01")))
      cat("\n")
      cat("  - Dataset end date:", format(as.Date(tail(dates,1), origin = "1970-01-01")))
      cat("\n")
      if (reinit){
        cat("  - Reinitialization date:",format(as.Date(reinit.date, origin = "1970-01-01")))
        cat("\n")
        cat("  - Use presample info:", use.presample.info)
        cat("\n")
      }
      cat("  - Model States and Standard Errors\n")
      base::print(out)
    },
    print = function() {
      reinit<-!is.null(reinit.date)
      out <- output(.self$estimate()) #KFS object
      if(is.null(q)){
        qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      }
      cat("SSModelDynamicGompertz Model")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("\n")
      cat("Cumulated Variable:\n")
      base::print(head(Y))
      cat("Number of observations:", length(.self$Y))
      cat("\n")
      cat("Signal-to-Noise Ratio (q):", 
          ifelse(is.null(q), paste(signif(qest,5), "(estimated)"), 
                 paste(q, ("(user specified)"))), "\n")
      cat("Seasonal components?",
          ifelse(is.null(seasonalComp(out)),
                 "No","Yes"),"\n")
      if (!is.null(reinit.date)){
        cat("Reinit date:",format(as.Date(reinit.date, origin = "1970-01-01")))
        cat("\n")
        cat("Use presample info:", use.presample.info)
      }
    },
    plot_diff =function(...){
      plot(diff(Y), ...)
    }
  )
)
