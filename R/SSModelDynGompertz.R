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
#' @title Class for Dynamic Gompertz Curve State-Space Model
#'
#' @description Class for Dynamic Gompertz Curve State-Space Model Object, which encapsulates
#' model settings and provides methods to obtain a FilterResults object and plot time series.
#' 
#' The dynamic Gompertz model with an integrated random walk (IRW) trend is defined as:
#' \deqn{\ln g_{t}= \delta_{t} + \varepsilon_{t}, \quad
#' \varepsilon_{t} \sim NID(0, \sigma_{\varepsilon}^{2}), \quad
#' t=2, ..., T,}
#' where \eqn{Y_t} is the cumulative variable, \eqn{y_t = \Delta Y_t}, and
#' \deqn{\ln g_{t} = \ln y_{t} - \ln Y_{t-1}.}
#' The trend component follows:
#' \deqn{\delta_{t} = \delta_{t-1} + \gamma_{t-1},}
#' \deqn{\gamma_{t} = \gamma_{t-1} + \zeta_{t}, \quad
#' \zeta_{t} \sim NID(0, \sigma_{\zeta}^{2}).}
#' Here, the observation disturbances \eqn{\varepsilon_{t}} and slope disturbances \eqn{\zeta_{t}} are independent and normally distributed. The signal-to-noise ratio,
#' \eqn{q_{\zeta} = \sigma_{\zeta}^{2} / \sigma_{\varepsilon}^{2}},
#' determines how rapidly the slope adjusts to new observations—higher values lead to faster changes, while lower values induce smoothness.
#' For models without seasonal terms (\code{sea.period = 0}), the priors are given by:
#' \deqn{\begin{pmatrix} \delta_1 \ \gamma_1 \end{pmatrix}
#' \sim N(a_1, P_1).}
#'
#' The diffuse prior is defined as \eqn{P_1 = \kappa I_{2\times 2}} with \eqn{\kappa \to \infty}, implemented via the \code{KFAS} package (Helske, 2017). For models with a seasonal component (\code{sea.period>1}), the prior mean vector \eqn{a_1} and prior covariance matrix \eqn{P_1} are extended accordingly.
#'
#' See the vignette for details on the variance matrix \eqn{Q} and the observation noise variance \eqn{Q}, \eqn{H = \sigma^2_{\varepsilon}}.
#' 
#' This class also supports the implementation of the reinitialisation
#' procedure, described in the vignette and also summarised below.
#' Let \eqn{t=r} denote the re-initialization date and \eqn{r_0} denote the
#' date at which the cumulative series is set to 0. As the growth rate of
#' cumulative cases is defined as \eqn{g_t\equiv \frac{y_t}{Y_{t-1}}}, we have:
#' \deqn{\ln g_t = \ln y_t - \ln Y_{t-1} \;\;\;\; t=1, \ldots, r}
#' \deqn{\ln g_t^r = \ln y_t - \ln Y_{t-1}^r \;\;\;\; t=r+1, \ldots, T}
#' \deqn{Y_{t}^{r}=Y_{t-1}^{r}+y_{t}  \;\;\;\; t=r,\ldots,T}
#' where \eqn{Y_{t}^{r}} is the cumulative cases after re-initialization. We
#' choose to set the cumulative cases to zero at \eqn{r_0=r-1, Y_{r-1}^{r}=0},
#' such that the growth rate of cumulative cases is available from \eqn{t=r+1}
#' onwards.
#' We reinitialise the model by specifying the prior distribution for the
#' initial states appropriately. See the vignette for details.
#' 
#' @field Y The cumulated variable.
#' @field q The signal-to-noise ratio (ratio of slope to irregular
#'   variance). Defaults to \code{'NULL'}, in which case no
#'   signal-to-noise ratio will be imposed. Instead, it will be estimated.
#'@field sea.period A positive integer specifying the period of seasonality used in the
#'   trigonometric seasonal component of the model. For example, use \code{7} for daily 
#'   data to model day-of-the-week effects. A value of \code{0} disables the seasonal 
#'   component entirely. The default is \code{7}, which is suitable for capturing 
#'   weekly seasonality in daily time series.
#'@field reinit.date (Only needed for reinitialization.) The reinitialisation date \eqn{r}. Should be
#' specified as an object of class \code{"Date"}. Defaults to NA, which 
#' represents the non-reinitialized version.
#' @field original.results (Only needed for reinitialization.) Rather than re-estimating the model up
#' to the \code{reinit.date}, a \code{FilterResults} class object can be
#' specified here and the parameters for the reinitialisation will be taken
#' from this object. Default is \code{NULL}. This parameter is optional.
#' @field use.presample.info (Only needed for reinitialization.) Logical value denoting whether or
#' not to use information from before the reinitialisation date in the
#' reinitialisation procedure. Default is \code{TRUE}. If \code{FALSE}, the
#' model is estimated from scratch from the reinitialisation date and no
#' attempt to use information from before the reinitialisation date is made.
#' @field xpred An \code{xts} object containing the dataset of exogenous variables 
#' to include in the model. Defaults to \code{NULL}.
#' @field ar1 Logical value indicating whether an ar1 component should be 
#' included in the model. Default is \code{FALSE}.
#' 
#' @importFrom xts periodicity last
#' @importFrom methods new
#' @importFrom xts xts
#' @importFrom zoo index
#' @importFrom KFAS SSModel fitSSM KFS
#' @importFrom magrittr %>%
#' @importFrom KFAS SSMtrend SSMseasonal
#' @import ggplot2
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' 
#' # Show summary of the model object
#' summary(model)
#' 
#' # Print a short description of the model object
#' print(model)
#' 
#' # Plot the time series in the model object
#' plot(model,title="Daily COVID cases in Gauteng", series.name="cases", MA=TRUE)
#' 
#' # Estimate a specified model
#' res <- estimate(model)
#' res
#'
#' @export SSModelDynamicGompertz
#' @exportClass SSModelDynamicGompertz
SSModelDynamicGompertz <- setRefClass(
  "SSModelDynamicGompertz",
  fields = list(
    Y = "xts",
    q = "ANY",
    sea.period="numeric",
    reinit.date = "ANY",
    original.results = "ANY",
    use.presample.info = "ANY",
    xpred="ANY",
    ar1="logical"
  ),
  methods = list(initialize = function(Y, q = NULL, 
                                       sea.period = 7,reinit.date=NULL, 
                                       original.results=NULL,
                                       use.presample.info=TRUE, xpred=NULL, 
                                       ar1=FALSE)
  {
    "Create an instance of the \\code{SSModelDynamicGompertz} class. Parameters 
    are defined in `fields` section. 
      \\subsection{Usage}{\\code{SSModelDynGompertzReinit$new(y, q = 0.005,
      reinit.date = as.Date(\"2021-05-12\",format = date.format))}}"
    if (!is.numeric(sea.period) || sea.period==1 || sea.period<0){
      stop("sea.period must be a positive integer that is not 1.")
    }
    Y <<- Y
    q <<- q
    sea.period <<- sea.period
    reinit.date <<- reinit.date
    original.results <<- original.results
    use.presample.info <<- use.presample.info
    xpred<<-xpred[index(Y)]
    ar1<<-ar1
  },
  estimate = function() {
    "Estimates the dynamic Gompertz curve model when applied to an object of
      class \\code{SSModelDynamicGompertz}.
      \\subsection{Return Value}{An object of class \\code{FilterResults}
      containing the result output for the estimated dynamic Gompertz curve
      model.}
      "
    update = function(pars, model, q) {
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
      }}
      \\subsection{Return Value}{\\code{KFS} model object.}"
      estH <- any(is.na(model$H))
      estQ <- any(is.na(model$Q))
      if ((!estH) & (!estQ)) {
        # If nothing to update then return model
        return(model)
      } else {
        nparQ <- 0
        # 1. Set seasonal noise
        if (estQ) {
          Q <- as.matrix(model$Q[, , 1])
          # Update diagonal elements
          naQd <- which(is.na(diag(Q)))
          if (ar1) {
            i.ar1 <- nrow(Q)
            naQd <- setdiff(naQd, i.ar1)
          }
          
          if (sea.period >1){
            nparQ <- 1
            Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
            diag(Q)[naQd] <- exp(0.5 * pars[nparQ])
            # Check for off-diagonal elements and raise error if found.
            naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
            if (length(naQnd) > 0) {
              stop("NotImplmentedError: Unexpected off-diagonal element updating")
            }
          }
        
          # 2. Set observation noise
          H <- as.matrix(model$H[, , 1])
          if (estH) {
            naHd <- which(is.na(diag(H)))
            H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
            nparQ<-nparQ+1
            diag(H)[naHd] <- exp(0.5 * pars[nparQ])
            model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
          }
            
        # 3. Set slope noise
        # Get index of slope, 1 before the seasonal component.
        model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
        i.slope <- 2
        # Estimate slope if no signal to noise ratio specified.
        if (is.null(q)) {
          nparQ<-nparQ+1
          Q.slope <- exp(0.5 * pars[nparQ])
        } else {
          Q.slope <- crossprod(H[naHd, naHd]) * q
        }
        model$Q[i.slope, i.slope, 1] <- Q.slope
        
        # 4. Set AR1 noise
        if (ar1){
          nparQ<-nparQ+1
          i.ar1 <- nrow(Q)
          Q[i.ar1, i.ar1] <- exp(0.5 * pars[nparQ])
          model$Q[i.ar1, i.ar1, 1] <- Q[i.ar1, i.ar1]
          
          nparQ<-nparQ+1
          T <- model$T[,,1]
          model$T[nrow(T),ncol(T),1] <- pars[nparQ]
        }
        }
      }
      return(model)
    }
    
    get_model = function(y,xpred=NULL){
      get_dynamic_gompertz_model = function(
    y,
    xpred,
    a1 = NULL,
    P1 = NULL,
    Q = NULL,
    H = NULL,
    T=NULL,
    R=NULL,
    newZ=NULL)
      { "Obtain the model object which is then used for 
        estimation."
        Ht <- if (is.null(H)) { NA } else { H }
        Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
        if (sea.period>1){
          Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
        }
        Qt.ar1 <- if (is.null(Q)) { NA } else {Q[dim(Q)[1],dim(Q)[2]]}
        
        # 1. Set prior on state as ~ N(a1, P1) if a1 supplied.
        use.prior <- if (!is.null(a1)) { TRUE } else { FALSE }
        
        # 2. Check whether there are exogenous predictors in model
        need.xpred<-!is.null(xpred)
        
        if (ar1){
          # 3. When needed, extract the AR1 coefficient
          ar1_coeff<-T[dim(T)[1],dim(T)[2]]
        }
        
        #Write out the model depending on case
        if (use.prior) {
          seasonal_idx<-grep("sea_trig", rownames(a1))
          trend_idx<-c(grep("level", rownames(a1)), 
                       grep("slope", rownames(a1)))
          #Case 1: With prior info, seasonality, xpred
          if (sea.period>1) {
            if (need.xpred){
              if (ar1){
                ss_model <-SSModel(
                  as.matrix(y) ~ SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope)),
                      a1 = a1[trend_idx],
                      P1 = P1[trend_idx, trend_idx]
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric",
                      a1 = a1[seasonal_idx],
                      P1 = P1[seasonal_idx, seasonal_idx]
                    )+SSMcustom(Z=1,T=ar1_coeff,R=1,Q=Qt.ar1, 
                                a1=a1[dim(a1)[1]], 
                                P1=P1[dim(a1)[1],dim(a1)[1]], 
                                state_names="ar1")
                  +SSMregression(~xpred),
                  H = Ht)
              } else {
                ss_model <-SSModel(
                  as.matrix(y) ~SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope)),
                      a1 = a1[trend_idx],
                      P1 = P1[trend_idx, trend_idx]
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric",
                      a1 = a1[seasonal_idx],
                      P1 = P1[seasonal_idx, seasonal_idx])
                  +SSMregression(~xpred),
                  H = Ht)
              }
            } else {
              #Case 2: With prior info, seasonality, no xpred
              if (ar1){
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope)),
                      a1 = a1[1:2],
                      P1 = P1[1:2, 1:2]
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric",
                      a1 = a1[3:(dim(a1)[1]-1)],
                      P1 = P1[3:(dim(a1)[1]-1), 3:(dim(a1)[1]-1)]
                    )+SSMcustom(Z=1,T=ar1_coeff,R=1,Q=Qt.ar1, 
                                a1=a1[dim(a1)[1]], 
                                P1=P1[dim(a1)[1],dim(a1)[1]], 
                                state_names="ar1"),
                  H = Ht)
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope)),
                      a1 = a1[1:2],
                      P1 = P1[1:2, 1:2]
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric",
                      a1 = a1[3:dim(a1)[1]],
                      P1 = P1[3:dim(a1)[1], 3:dim(a1)[1]]),
                  H = Ht
                ) 
              }
            }
          } else {
            #Case 3: With prior info, no seasonality, yes xpred
            if (need.xpred){
              if (ar1){
                 ss_model <-SSModel(
                  as.matrix(y) ~ SSMtrend(
                    degree = 2,
                    Q = list(matrix(0), matrix(Qt.slope)),
                    a1 = a1[trend_idx],
                    P1 = P1[trend_idx, trend_idx]
                  )+SSMcustom(Z=1,T=ar1_coeff,R=1,Q=Qt.ar1, 
                                a1=a1[dim(a1)[1]], 
                                P1=P1[dim(a1)[1],dim(a1)[1]], 
                                state_names="ar1")
                  +SSMregression(~xpred),
                  H = Ht)
              } else {
                ss_model <-SSModel(
                  as.matrix(y) ~ SSMtrend(
                    degree = 2,
                    Q = list(matrix(0), matrix(Qt.slope)),
                    a1 = a1[trend_idx],
                    P1 = P1[trend_idx, trend_idx])
                  +SSMregression(~xpred),
                  H = Ht)
              }
            } else {
              #Case 4: With prior info, no seasonality, no xpred
              if (ar1){
                ss_model <- SSModel(
                  as.matrix(y) ~SSMtrend(
                    degree = 2,
                    Q = list(matrix(0), matrix(Qt.slope)),
                    a1 = a1[1:2],
                    P1 = P1[1:2, 1:2])
                  +SSMcustom(Z=1,T=ar1_coeff,R=1,Q=Qt.ar1, 
                             state_names="ar1"),
                  H = Ht)
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~SSMtrend(
                    degree = 2,
                    Q = list(matrix(0), matrix(Qt.slope)),
                    a1 = a1[1:2],
                    P1 = P1[1:2, 1:2]),
                  H = Ht)
              }
              }
          } 
          n.pars <- 0
        } else {
          #Case 5: No prior info, yes seasonality, yes xpred
          if (need.xpred){
            if (sea.period>1) {
              if(ar1){
                ss_model <- SSModel(
                  as.matrix(y) ~SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric")+
                    SSMcustom(Z=1,T=1,R=1,Q=matrix(NA),state_names="ar1")
                  +SSMregression(~xpred),
                  H = matrix(Ht)
                )
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric")
                  +SSMregression(~xpred),
                  H = matrix(Ht)
                )
              }
          #Case 6: No prior info, no seasonality, yes xpred
            } else {
              if (ar1){
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    )+SSMcustom(Z=1,T=1,R=1,Q=matrix(NA),state_names="ar1")
                  +SSMregression(~xpred),
                  H = matrix(Ht)
                )
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    )+SSMregression(~xpred),
                  H = matrix(Ht)
                )
              }
            } 
          } else {
            #Case 7: No prior info, yes seasonality, no xpred
            if (sea.period>1) {
              if (ar1){
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric")+
                    SSMcustom(Z=1,T=1,R=1,Q=matrix(NA),state_names="ar1"),
                  H = matrix(Ht)
                )
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    ) +
                    SSMseasonal(
                      period = sea.period,
                      Q = Qt.seas,
                      sea.type = "trigonometric"),
                  H = matrix(Ht)
                )
              }
            } else {
              #Case 8: No prior info, no seasonality, no xpred
              if (ar1){
                ss_model <- SSModel(
                  as.matrix(y) ~SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    )+
                    SSMcustom(Z=1,T=1,R=1,Q=matrix(NA),state_names="ar1"),
                  H = matrix(Ht)
                )
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))),
                  H = matrix(Ht))
              }
            } 
          }
          n.pars <- sum(is.na(ss_model$Q)) + sum(is.na(ss_model$H))
          if (!is.null(q)){n.pars<-n.pars-1}
        }
        if (ar1){
          out <- list(model = ss_model, inits = c(rep(0,n.pars),1))
        } else {
          out <- list(model = ss_model, inits = rep(0, n.pars))
        }
        return(out)
      }
      
      if (is.null(reinit.date)){
        model <- get_dynamic_gompertz_model(
          y, xpred=xpred
        )
        return(model)
      } else{
        #Select relevant xpred
        xpred1<-xpred[zoo::index(xpred) <= reinit.date]
        xpred2<-xpred[zoo::index(xpred) > reinit.date]
        
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
            model <- SSModelDynamicGompertz$new(Y = Y[idx.est],
                                                sea.period=sea.period, 
                                                xpred=xpred1, q = q, ar1=ar1)
            res.original <- model$estimate()
            model_output <- output(res.original)
          } else {
            model_output <- output(original.results)
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
          a1 <- NULL; P1 <- NULL; Qt <- NULL; Ht <- NULL; Tt<- NULL
        }
        out <- get_dynamic_gompertz_model(
          y = y.reinit, xpred=xpred2,
          a1 = a1, P1 = P1, Q = Qt, H = Ht, T=Tt)
        
        out[['index']] <- index(y.reinit)
        return(out)
      }
    }
    
    # 1. Get LDL of cumulative series Y.
    y <- tsgc::df2ldl(Y)
    
    # 2. Add update / model methods
    updatefn <- purrr::partial(
      update, ... =, q = q
    )
    
    model <- get_model(y, xpred=xpred)

    # 3. Estimate via MLE unknown params
    model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                        method = 'BFGS')
    
    # 4. Run smoother/filter
    model_output <- KFS(model_fit$model)
    
    # 5. Get truncated index from model if using a reinitialisation in model
    date.index <- if (!is.null(model$index)) { model$index } else { index(y) }
    
    results <- FilterResults$new(
      data_xts = Y,
      need.xpred=!is.null(xpred),
      index = date.index,
      reinit.date=reinit.date,
      ar1=ar1,
      sea.period=sea.period,
      output = model_output
    )
    return(results)
  },
  summary = function() {
    "Supplies details of the SSModelDynamicGompertz object, such as estimated 
      parameter values, start and end dates of estimation."
    out <- suppressWarnings(output(.self$estimate()))
    if(is.null(q)){
      qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
    }
    reinit<-!is.null(reinit.date)
    dates<-index(Y)
    ar1_comp<-matrixKFS(out,"T")["ar1","ar1",1]
    
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
    if (ar1){
      cat("AR(1) coefficient:", signif(ar1_comp,3))
      cat("\n")
    }
    cat("Model Details:\n")
    cat("  - Model Type: Dynamic Gompertz Curve")
    if (reinit) {
      cat(" (Reinitialized)")
    }
    cat("\n")
    cat("  - Seasonal Component: ", ifelse(sea.period>1, "Trigonometric", "None"), "\n")
    cat("  - Period of Seasonality: ", ifelse(sea.period>1, sea.period, "N/A"), "\n")
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
    if (!is.null(xpred)){
      cat("  - Exogenous predictors dataset")
      base::print(head(.self$xpred))
    }
    cat("  - Model States and Standard Errors\n")
    base::print(out)
  },
  print = function() {
    "Provides a quick description of SSModelDynamicGompertz object, providing 
      model states and standard errors."
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
    cat("Exogenous predictors?", ifelse(is.null(xpred),
                                         "No","Yes"),"\n")
    if (!is.null(reinit.date)){
      cat("Reinit date:",format(as.Date(reinit.date, origin = "1970-01-01")))
      cat("\n")
      cat("Use presample info:", use.presample.info)
    }
  },
  plot =function(title=NULL, series.name="target variable", MA=TRUE){
    "Plots the lagged differences of the cumulated dataset \\code{Y} in this 
      \\code{SSModelDynamicGompertz} object against time, which could represent 
      daily cases.
      \\subsection{Parameters}{\\itemize{
        \\item{\\code{MA} A logical value indicating whether 7-day centered moving 
        average should be plotted. Defaults to \\code{TRUE}.
        }
        \\item{\\code{title} Title for forecast plot. Enter as text string. 
        \\code{NULL} (i.e. no title) by default.}
        \\item{\\code{series.name} The name of the series the growth rate is being computed
        for. E.g. \\code{'cases'}. Default is `target variable`.}
}
      }
      \\subsection{Return Value}{A plot of the lagged differences of the 
      cumulated dataset \\code{Y} against time.}
      "
  cumulative_cases <- Y  
  
  # Calculate a centred 7-day moving average of daily differences.
  ma.cent.new.cases <- zoo::rollmean(diff(cumulative_cases), 7, align = "center")
  
  # Identify the date with maximum new cases.
  ma.cent.wave.3.idx.max <- tsgc::argmax(ma.cent.new.cases) %>% zoo::index()
  
  # Prepare data for plotting by combining actual new cases and the moving average.
  d <- cbind(diff(cumulative_cases), ma.cent.new.cases)
  colnames(d) <- c('New Cases', 'Centered 7-day MA')
  d.df <- data.frame(
    Date = index(d),
    New.Cases = coredata(d[, 1]),
    Centered.7.day.MA = coredata(d[, 2])
  )
  
  # Create base plot
  data_plot <- ggplot(data = d.df, aes(x = Date)) +
    geom_line(aes(y = New.Cases, color = "New Cases"), linewidth = 0.1) +
    scale_y_continuous(n.breaks = 10) +
    labs(x = "Date", y = paste("New", series.name), title = title)+
    scale_x_date(date_breaks = "60 days") +
    theme_light(base_size = 12) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.2, 0.85),
      legend.title = element_text(size = 2),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold")
    )
  
  # Conditionally add Centered 7-day MA line
  if (MA) {
    data_plot <- data_plot + 
      geom_line(aes(y = Centered.7.day.MA, color = "Centered 7 day MA"), linewidth = 1)+ 
      scale_color_manual(
        name = '',
        values = c('Centered 7 day MA' = 'red')
      )
  } 
  data_plot
  }
  )
)


