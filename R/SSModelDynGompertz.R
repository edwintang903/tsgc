setOldClass("KFS")
setOldClass("idx_series")

#' @title Class for designing a Dynamic Gompertz Curve State-Space Model
#'
#' @description Class for Dynamic Gompertz Curve State-Space Model Object, which encapsulates
#' model settings and provides methods to obtain a FilterResults object.
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
#' See the vignette for details on the state disturbance variance matrix \eqn{Q} and the observation noise variance \eqn{H = \sigma^2_{\varepsilon}}.
#' 
#' This class also supports the implementation of the reinitialisation
#' procedure, described in the vignette and also summarised below.
#' Let \eqn{t=r} denote the re-initialization position and \eqn{r_0} denote the
#' position at which the cumulative series is set to 0. As the growth rate of
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
#' @field Y The cumulated variable, as an \code{idx_series}. Must be
#'   strictly increasing.
#' @field q The signal-to-noise ratio (ratio of slope to irregular
#'   variance). Defaults to \code{'NULL'}, in which case no
#'   signal-to-noise ratio will be imposed. Instead, it will be estimated.
#' @field sea.period A positive integer specifying the period of seasonality used in the
#'   trigonometric seasonal component of the model. For example, use \code{7} for daily 
#'   data to model day-of-the-week effects. A value of \code{0} disables the seasonal 
#'   component entirely. The default is \code{7}, which is suitable for capturing 
#'   weekly seasonality in daily time series.
#' @field reinit.idx (Only needed for reinitialization.) The
#' reinitialisation position \eqn{r}, as a single integer. Defaults to
#' \code{NULL}, which represents the non-reinitialized version.
#' @field original.results (Only needed for reinitialization.) Rather than re-estimating the model up
#' to the \code{reinit.idx}, a \code{FilterResults} class object can be
#' specified here and the parameters for the reinitialisation will be taken
#' from this object. Default is \code{NULL}. This parameter is optional.
#' @field use.presample.info (Only needed for reinitialization.) Logical value denoting whether or
#' not to use information from before the reinitialisation position in the
#' reinitialisation procedure. Default is \code{TRUE}. If \code{FALSE}, the
#' model is estimated from scratch from the reinitialisation position and no
#' attempt to use information from before the reinitialisation position is made.
#' @field xpred An \code{idx_series} object containing the dataset of exogenous
#' variables to include in the model. Defaults to \code{NULL}.
#' @field ar1 Logical value indicating whether an ar1 component should be 
#' included in the model. Default is \code{FALSE}.
#' @field start Integer position marking the start of the estimation period.
#' @field end Integer position marking the end of the estimation period.
#' @field calendar An optional \code{idx_calendar} object describing how the
#' integer positions used by \code{Y} (and the rest of this model's
#' \code{idx_series} fields) map back to calendar time. Purely cosmetic:
#' never consulted by estimation/filtering, only carried along so that
#' plotting can later translate positions to dates. Defaults to
#' \code{NULL}, in which case plots will be labelled with raw integer
#' positions instead of dates.
#' 
#' @importFrom methods new setRefClass setOldClass
#' @importFrom KFAS SSModel fitSSM KFS SSMtrend SSMseasonal SSMregression SSMcustom
#' @importFrom magrittr %>%
#'
#' @examples
#' library(tsgc)
#' set.seed(1)
#' Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1)
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 100)
#' 
#' # Show summary of the model object
#' summary(model)
#' 
#' # Print a short description of the model object
#' print(model)
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
    Y = "idx_series",
    q = "ANY",
    sea.period="numeric",
    reinit.idx = "ANY",
    original.results = "ANY",
    use.presample.info = "ANY",
    xpred="ANY",
    ar1="logical",
    start="ANY",
    end="ANY",
    calendar="ANY"),
  methods = list(initialize = function(Y, q = NULL, 
                                       sea.period = 7,reinit.idx=NULL, 
                                       original.results=NULL,
                                       use.presample.info=TRUE, xpred=NULL, 
                                       ar1=FALSE, start=idx_range(Y)[1], 
                                       end=idx_range(Y)[2], calendar=NULL)
  {
    "Create an instance of the \\code{SSModelDynamicGompertz} class. Parameters 
    are defined in `fields` section. 
      \\subsection{Usage}{\\code{SSModelDynamicGompertz$new(Y = y, q = 0.005,
      reinit.idx = 45)}}"
    if (length(sea.period) != 1 || 
        !isTRUE(all.equal(sea.period, as.integer(sea.period)))||
        sea.period==1 || sea.period<0){
      stop("sea.period must be a non-negative integer that is not 1.")
    } 
    if (!is.null(original.results) && !inherits(original.results, "FilterResults")){
      stop("original.results must be NULL or an object of class FilterResults.")
    }
    if (!is.null(xpred) && !is_idx_series(xpred)){
      stop("xpred must be NULL or an idx_series object.")
    } 
    if (!is.null(calendar) && !is_idx_calendar(calendar)){
      stop("calendar must be NULL or an idx_calendar object.")
    }
    Y <<- get_timeframe(Y,start,end)
    q <<- q
    sea.period <<- sea.period
    reinit.idx <<- reinit.idx
    original.results <<- original.results
    use.presample.info <<- use.presample.info
    xpred<<-get_timeframe(xpred,start,end)
    ar1<<-ar1
    start<<-start
    end<<-end
    calendar<<-calendar
  },
  estimate = function() {
    "Estimates the dynamic Gompertz curve model when applied to an object of
      class \\code{SSModelDynamicGompertz}.
      \\subsection{Return Value}{An object of class \\code{FilterResults}
      containing the result output for the estimated dynamic Gompertz curve
      model.}
      "
    if (any(na.omit(idx_values(idx_diff(Y, 1L)))<=0)){
      stop("Y must be strictly increasing. If the cumulative 
           values exhibit plateaus it is necessary to add small increments to 
           eliminate flat segments and allow model estimation. This can be done 
           by ensuring the non-cumulated series is strictly positive.")
    }
    
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
            model$Q[i.slope, i.slope, 1] <- crossprod(Q.slope)
          } else {
            model$Q[i.slope, i.slope, 1] <- crossprod(H[naHd, naHd]) * q
          }
          
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
    xreg,
    a1 = NULL,
    P1 = NULL,
    Q = NULL,
    H = NULL,
    T=NULL,
    R=NULL,
    newZ=NULL)
      { "Obtain the model object which is then used for 
        estimation."
        # Named `xreg` rather than `xpred` to avoid shadowing the
        # SSModelDynamicGompertz RefClass field `xpred`. KFAS's
        # SSMregression(~xreg) names the fitted coefficient state after
        # this formula variable name, so it is also user-visible in
        # print()/summary() output.
        if (is_idx_series(xreg) && is_idx_series(y)) {
          xreg <- xreg[idx_positions(y)]
        }
        if (is_idx_series(xreg)) { xreg <- idx_values(xreg) }
        if (is_idx_series(y)) { y <- as.matrix(idx_values(y)) }
        Ht <- if (is.null(H)) { NA } else { H }
        Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
        if (sea.period>1){
          Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
        }
        Qt.ar1 <- if (is.null(Q)) { NA } else {Q[dim(Q)[1],dim(Q)[2]]}
        
        # 1. Set prior on state as ~ N(a1, P1) if a1 supplied.
        use.prior <- if (!is.null(a1)) { TRUE } else { FALSE }
        
        # 2. Check whether there are exogenous predictors in model
        need.xpred<-!is.null(xreg)
        
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
                  +SSMregression(~xreg),
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
                  +SSMregression(~xreg),
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
                  +SSMregression(~xreg),
                  H = Ht)
              } else {
                ss_model <-SSModel(
                  as.matrix(y) ~ SSMtrend(
                    degree = 2,
                    Q = list(matrix(0), matrix(Qt.slope)),
                    a1 = a1[trend_idx],
                    P1 = P1[trend_idx, trend_idx])
                  +SSMregression(~xreg),
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
                  +SSMregression(~xreg),
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
                  +SSMregression(~xreg),
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
                  +SSMregression(~xreg),
                  H = matrix(Ht)
                )
              } else {
                ss_model <- SSModel(
                  as.matrix(y) ~
                    SSMtrend(
                      degree = 2,
                      Q = list(matrix(0), matrix(Qt.slope))
                    )+SSMregression(~xreg),
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
      
      if (is.null(reinit.idx)){
        model <- get_dynamic_gompertz_model(
          y, xreg=xpred
        )
        return(model)
      } else{
        if (is.null(xpred)) {
          xpred1 <- NULL
          xpred2 <- NULL
        } else {
          #Select relevant xpred
          xpred_pos <- idx_positions(xpred)
          xpred1 <- if (any(xpred_pos <= reinit.idx)) xpred[xpred_pos[xpred_pos <= reinit.idx]] else NULL
          xpred2 <- if (any(xpred_pos > reinit.idx)) xpred[xpred_pos[xpred_pos > reinit.idx]] else NULL
        }
        
        # 4.1. Position for reinitialisation, t_0
        stopifnot(reinit.idx %in% idx_positions(Y))
        # Coerce to a plain numeric scalar: idx_values() on a single-row
        # idx_series can return a 1x1 matrix, which does not broadcast
        # against the n x 1 matrix idx_values(lag.Y) below.
        Y.t.r_0 <- as.numeric(idx_values(Y[reinit.idx - 1]))
        
        # 4.2 Reinitialisation:
        #   ln g_t^r = ln g_t + ln (Y_{t-1}/Y_{t-1}^r), where Y_t^r=Y_t-Y_{r_0}.
        y_pos <- idx_positions(y)
        reinit_pos <- y_pos[y_pos > reinit.idx]
        lag.Y <- idx_lag(Y, 1L)[reinit_pos]
        # Coerce to plain numeric vectors before arithmetic, since y$data
        # and Y$data may differ in shape (plain vector vs matrix).
        y.vals <- as.numeric(idx_values(y[reinit_pos]))
        lag.Y.vals <- as.numeric(idx_values(lag.Y))
        y.reinit <- idx_series(
          y.vals + log(lag.Y.vals / (lag.Y.vals - Y.t.r_0)),
          start = reinit_pos[1]
        )
        
        # 4.3 Run Kalman filter/smoother on new series with non-diffuse prior
        if (use.presample.info) {
          # Either estimate full model here or take results from previous model.
          if (is.null(original.results)) {
            # NB. Restrict sample to t<=r - position of reinitialisation.
            model <- SSModelDynamicGompertz$new(Y = Y,
                                                sea.period=sea.period, 
                                                xpred=xpred1, q = q, ar1=ar1,
                                                start=start,
                                                end=reinit.idx,
                                                calendar=calendar)
            res.original <- model$estimate()
            model_output <- output(res.original)
          } else {
            model_output <- output(original.results)
          }
          
          # 4.3 Reset slope to 0 and add constant to initial value for level.
          # where reinit.idx is t=r
          idx <- which(reinit.idx == idx_positions(y))
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
          a1["level",] <- a1["level",] + log(idx_values(Y[reinit.idx]) / (idx_values(Y[reinit.idx]) - Y.t.r_0))
          
        } else {
          # Don't use presample info
          a1 <- NULL; P1 <- NULL; Qt <- NULL; Ht <- NULL; Tt<- NULL
        }
        out <- get_dynamic_gompertz_model(
          y = y.reinit, xreg=xpred2,
          a1 = a1, P1 = P1, Q = Qt, H = Ht, T=Tt)
        
        out[['index']] <- idx_positions(y.reinit)
        return(out)
      }
    }
    
    # 1. Get LDL of cumulative series Y.
    y <- df2ldl(Y)
    
    # 2. Obtain the SSModel 
    model <- get_model(y, xpred=xpred)
    
    # 3. Add update methods to enforce signal-to-noise ratio
    updatefn <- purrr::partial(update, ... =, q = q)
    
    # Estimate via MLE unknown params
    model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                        method = 'BFGS')
    
    # 4. Run smoother/filter
    model_output <- KFS(model_fit$model)
    
    # 5. Get truncated index from model if using a reinitialisation in model
    idx.positions <- if (!is.null(model$index)) { model$index } else { idx_positions(y) }
    
    results <- FilterResults$new(
      data = Y,
      xpred_logical = !is.null(xpred),
      index = idx.positions,
      reinit.idx =reinit.idx,
      ar1=ar1,
      sea.period=sea.period,
      output = model_output,
      calendar = calendar
    )
    return(results)
  },
  summary = function() {
    "Estimates the model and delegates to the resulting FilterResults
      object's summary() method, so that summary(model) reports the same
      fitted-results summary as summary(res)."
    .self$estimate()$summary()
  },
  print = function() {
    "Provides a quick description of the SSModelDynamicGompertz object, providing 
      model states and standard errors."
    reinit<-!is.null(reinit.idx)
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
    base::print(head(idx_values(Y)))
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
    if (!is.null(reinit.idx)){
      cat("Reinit position:",reinit.idx)
      cat("\n")
      cat("Use presample info:", use.presample.info)
    }
  }
  )
)