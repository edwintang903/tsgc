# Created by: Craig Thamotheram
# Refactored: analysis-only FilterResults class operating on idx_series
# (integer-indexed) data. All plotting methods have been removed from this
# file; they will be reintroduced elsewhere as a purely cosmetic layer that
# translates integer positions back to calendar time before plotting.

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
setOldClass("idx_series")
#'
#' @title Class for estimated Dynamic Gompertz Curve model
#'
#' @description Class for estimated Dynamic Gompertz Curve model and contains
#' methods to extract smoothed/filtered estimates of the states, the level of
#' the incidence variable \eqn{y}, and forecasts of \eqn{y}. The output from the estimate method
#' of the SSModelDynamicGompertz class is of the class FilterResults.
#' 
#' @field data An \code{idx_series} object containing the non-reinitialized
#' cumulated variable.
#' @field xpred_logical Logical value indicating whether exogenous predictors were 
#' used to estimate the FilterResults object. 
#' @field index The integer positions of the observations used in
#' estimation (of \code{data}).
#' @field reinit.idx The reinitialisation position (a single integer) of the
#' estimated \code{SSModelDynamicGompertz} model (if applicable).
#' @field ar1 Logical value indicating whether an ar1 component should be 
#' included in the model.
#' @field output A \code{KFS} results object obtained after fitting a 
#' \code{SSModelDynamicGompertz} model.
#' @field xpred.new An \code{idx_series} object containing exogenous
#' predictors to be used in prediction. Defaults to \code{NULL}, and should
#' be provided if xpred is used for model estimation.
#' @field sea.period The period of seasonality, inherited from the estimated 
#' \code{SSModelDynamicGompertz} model. For a day-of-the-week
#'   effect with daily data, this would be 7. 
#' 
#' @references Harvey, A. C. and Kattuman, P. (2021). A Farewell to R:
#' Time Series Models for Tracking and
#' Forecasting Epidemics, Journal of the Royal Society Interface, vol 18(182):
#' 20210179
#'
#' @importFrom stats predict
#' @importFrom magrittr %>%
#' @importFrom methods new setRefClass setOldClass
#' @importFrom abind abind
#' 
#' @examples
#' library(tsgc)
#' set.seed(1)
#' Y <- idx_series(cumsum(rpois(120, 8)) + 1, start = 1)
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = Y, q = 0.005, end = 100)
#' 
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Show summary of object
#' summary(res)
#' 
#' # Print a short description of the object
#' print(res)
#' 
#' # Print estimation results
#' res$print_estimation_results()
#' 
#' # Forecast 7 periods ahead from the end of the estimation window
#' res$predict_level(n.ahead = 7,
#'   confidence.level = 0.68, sea.on=TRUE)
#'   
#' # Forecast 7 periods ahead from the model and return filtered states
#' res$predict_all(n.ahead = 7, return.all = TRUE)
#' 
#' # Return the filtered growth rate and its components
#' res$get_growth_y(return.components = TRUE)
#' 
#' # Return smoothed growth rate of incidence variable and its confidence
#' # interval
#' res$get_gy_ci(smoothed = TRUE, confidence.level = 0.68)
#'
#' # Return MAPE of forecast
#' res$mapes(n.ahead=7,Y)
#'
#' @export
#'
FilterResults <- setRefClass(
  "FilterResults",
  fields = list(
    data = "idx_series",
    xpred_logical = "ANY",
    xpred.new="ANY",
    index = "ANY",
    reinit.idx= "ANY",
    ar1 = "logical",
    output = "KFS",
    sea.period="numeric"),
  methods = list(
    initialize = function(data,xpred_logical,index,reinit.idx, ar1, 
                          output, sea.period, xpred.new=NULL)
    {
      "Create an instance of the \\code{FilterResults} class with fields defined
      earlier in the fields section."
      data<<-data
      index <<- index
      xpred_logical<<-xpred_logical
      xpred.new<<-xpred.new
      reinit.idx<<-reinit.idx
      ar1<<-ar1
      output <<- output
      sea.period<<-sea.period
    },
    predict_level = function(
    n.ahead,
    confidence.level=0.68,
    sea.on = TRUE, 
    return.diff=TRUE)
    {
      "Forecast the cumulated variable or the incidence of it. This function returns
      the forecast of the cumulated variable \\eqn{Y}, or the forecast of the incidence of the cumulated variable, \\eqn{y}. For
      example, in the case of an epidemic, \\eqn{y} might be daily new cases of
      the disease and
       \\eqn{Y} the cumulative number of recorded infections.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of periods ahead you wish to forecast from
        the end of the estimation window.}
        \\item{\\code{confidence.level} The confidence level for the log growth
         rate that should be used to compute
        the forecast intervals of \\eqn{y}.}
        \\item{\\code{sea.on} Logical value indicating whether to return the prediction 
        of just the trend or prediction incorporating seasonality. Deafults to \\code{TRUE}.}
        \\item{\\code{return.diff} Logical value indicating whether to return the cumulated variable,
        \\eqn{Y}, or the incidence of it,
        \\eqn{y} (i.e., the first difference of the cumulated variable). Default is
        \\code{TRUE}.}
      }}
      \\subsection{Return Value}{\\code{idx_series} object containing the point
      forecasts and upper and lower bounds of
      the forecast interval.}"
      if (!is.null(reinit.idx)){
        y.cum<-reinitialise_dataframe(data, reinit.idx)
      } else {
        y.cum<-data
      }
      model <- modelKFS(output)
      n <- attr(model, "n")
      p <- attr(model, "p")
      
      filtered.out <- .self$predict_all(n.ahead, sea.on = sea.on,
                                        return.all = FALSE, 
                                        confidence.level = confidence.level)
      
      # # 1. Extract parameters.
      timespan <- n + 0:n.ahead
      
      # Calculate g.t as exponent of y.t
      yhat_mat <- idx_values(gety.hat(filtered.out))
      g.t <- exp(yhat_mat[,1])
      g.t.lwr <- exp(yhat_mat[,2])
      g.t.upr <- exp(yhat_mat[,3])
      
      # Forecast positions: last position of estimation window through
      # n.ahead positions beyond it.
      last.pos <- idx_range(y.cum)[2]
      fc.positions <- seq.int(last.pos, length.out = n.ahead + 1)
      
      y.hat <- matrix(NA_real_, nrow = n.ahead + 1, ncol = 3)
      y.hat[1, 1] <- idx_values(y.cum[last.pos])
      for (i in seq_len(n.ahead)) {
        # Update level
        y.hat[i + 1, 1] <- y.hat[i, 1] * (1 + g.t[i])
        
        # Make prediction intervals
        y.hat[i + 1, 2] <- y.hat[i, 1] * (1 + g.t.lwr[i])
        y.hat[i + 1, 3] <- y.hat[i, 1] * (1 + g.t.upr[i])
      }
      y.hat <- idx_series(y.hat, start = fc.positions[1])
      
      # Difference output if requested
      d <- if (return.diff) { idx_diff(idx_series(y.hat$data[,1], start=y.hat$start), 1L) } else {
        idx_series(y.hat$data[-1, 1], start = y.hat$start + 1L)
      }
      
      ci_bounds <- if (return.diff) {
        base_mat <- y.hat$data[-1, 2:3, drop = FALSE] - y.hat$data[-nrow(y.hat$data), 1]
        base_mat + idx_values(d)
      } else {
        y.hat$data[-1, 2:3, drop = FALSE]
      }
      
      out <- idx_series(cbind(fit = idx_values(d), lower = ci_bounds[,1], upper = ci_bounds[,2]),
                        start = d$start)
      return(out)
    },
    print_estimation_results = function() {
      "Prints a table of estimated parameters in a format ready to paste into
      LaTeX."
      H <- output$model$H[, , 1]
      Q_gamma <- output$model$Q[2, 2, 1]
      has_seasonal <- sea.period > 1
      
      if (has_seasonal) {
        Q_seasonal <- output$model$Q[3, 3, 1]
        tbl <- data.frame(
          a = format(H, digits = 3),
          b = format(Q_gamma, digits = 3),
          c = format(Q_seasonal, digits = 3),
          d = format(Q_gamma / H, digits = 4))
        header.names <- c('$\\sigma_\\varepsilon^2$',
                          '$\\sigma_\\gamma^2$',
                          '$\\sigma_{seas}^2$',
                          'q')
      } else {
        tbl <- data.frame(
          a = format(H, digits = 3),
          b = format(Q_gamma, digits = 3),
          d = format(Q_gamma / H, digits = 4))
        header.names <- c('$\\sigma_\\varepsilon^2$',
                          '$\\sigma_\\gamma^2$',
                          'q')
      }
      
      out <- tbl %>%
        kableExtra::kbl(
          caption = "Estimated parameters",
          col.names = header.names,
          format = 'latex',
          booktabs = TRUE,
          escape = FALSE
        ) %>%
        kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria") %>%
        kableExtra::footnote(general = " ")
      
      return(out)
    },
    predict_all = function(n.ahead, sea.on = TRUE, return.all = FALSE, confidence.level = 0.68) {
      "Returns forecasts of the incidence variable \\eqn{y}, the state variables
       and the conditional covariance matrix
      for the states.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of forecasts you wish to create from
        the end of your sample period.}
        \\item{\\code{sea.on} Logical value indicating whether seasonal
        components should be included in the state-space model or not. Default is \\code{TRUE}.}
        \\item{\\code{return.all} Logical value indicating whether to return
        all filtered estimates and forecasts.
        (\\code{TRUE}) or only the forecasts (\\code{FALSE}). Default is
        \\code{FALSE}.}
        \\item{\\code{confidence.level} The confidence level for the log growth
         rate that should be used to compute. Confidence intervals only reported
         for the incidence variable \\eqn{y}.}
      }}
      \\subsection{Return Value}{\\code{idx_series} objects containing the forecast
      (and filtered, where applicable) level
      of \\eqn{y} (\\code{y.hat}), \\eqn{\\delta} (\\code{level.t.t}),
      \\eqn{\\gamma} (\\code{slope.t.t}), vector of states including the
      seasonals where applicable (\\code{a.t.t}) and covariance matrix of all
      states including seasonals where applicable (\\code{P.t.t}).}"
      
      new.model <- modelKFS(output)
      oldn<-attr(new.model, 'n')
      new.model$y <- rbind(
        gety(new.model),
        matrix(NA, ncol = ncol(gety(new.model)), nrow = n.ahead)) %>% as.ts()
      
      attr(new.model, 'n') <- as.integer(oldn + n.ahead)
      
      if (xpred_logical){
        if (is.null(xpred.new)){
          stop("xpred.new cannot be NULL.")
        } else {
          firstpred <- tail(index,1) + 1L
          
          xpred.new<<-get_timeframe(xpred.new, firstpred, firstpred + n.ahead - 1L)
          
          newZ<-array(new.model$Z[,,dim(new.model$Z)[3]], 
                      dim = c(dim(new.model$Z)[1], dim(new.model$Z)[2], n.ahead))
          xpred.new.mat <- as.matrix(idx_values(xpred.new))
          newZ[,1:dim(xpred.new.mat)[2],]<-t(xpred.new.mat)
          
          new.model$Z <- abind::abind(
            new.model$Z,
            newZ,
            along = 3
          )
          
          model_output <- KFS(new.model)
          new.Q <- new.model$Q
          xpred.new.vec <- xpred.new.mat
          if (ar1){
            #AR1 with sea.period
            if (sea.period > 1){
              ar1_index<-dim(new.Q)[1]
              newdata<-SSModel(rep(NA,dim(xpred.new.vec)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMseasonal(
                                 period = sea.period, 
                                 Q = new.Q[3,3,1],
                                 sea.type = "trigonometric")
                               +SSMregression(~xpred.new.vec)
                               +SSMcustom(Z=1,T=1,R=1,Q=new.Q[ar1_index,ar1_index,1],state_names="ar1"))
            } else {
              #AR1 and no sea.period
              ar1_index<-dim(new.Q)[1]
              newdata<-SSModel(rep(NA,dim(xpred.new.vec)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMregression(~xpred.new.vec)
                               +SSMcustom(Z=1,T=1,R=1,Q=new.Q[ar1_index,ar1_index,1],state_names="ar1"))
            }
          } else {
            #sea period only
            if (sea.period > 1){
              newdata<-SSModel(rep(NA,dim(xpred.new.vec)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMseasonal(
                                 period = sea.period, 
                                 Q = new.Q[3,3,1],
                                 sea.type = "trigonometric")
                               +SSMregression(~xpred.new.vec))
            } else {
              #no sea period 
              newdata<-SSModel(rep(NA,dim(xpred.new.vec)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMregression(~xpred.new.vec))
            }
          }
          
          if (sea.on == TRUE) {
            y.hat.kfas <- predict(
              output$model, interval = 'confidence', level = confidence.level,
              newdata = newdata, states = 'all')
            y.t.t <- predict(output$model, interval = 'confidence', 
                             level = confidence.level,
                             states = 'all')
          } else {
            y.hat.kfas <- predict(
              output$model, interval = 'confidence', level = confidence.level,
              newdata = newdata, states = c("level","regression","custom"))
            y.t.t <- predict(output$model, interval = 'confidence', 
                             level = confidence.level,
                             states = c("level","regression","custom"))
          }
        }
      } else {
        model_output <- KFS(new.model)
        
        if (sea.on == TRUE) {
          y.hat.kfas <- predict(
            output$model, interval = 'confidence', level = confidence.level,
            n.ahead = n.ahead, states = 'all')
          y.t.t <- predict(output$model, interval = 'confidence', 
                           level = confidence.level,
                           states = 'all')
        } else {
          y.hat.kfas <- predict(
            output$model, interval = 'confidence', level = confidence.level,
            n.ahead = n.ahead, states = c("level","regression","custom"))
          y.t.t <- predict(output$model, interval = 'confidence', 
                           level = confidence.level,
                           states = c("level","regression","custom"))
        }
      } 
      
      positions <- seq.int(index[1], length.out = (oldn + n.ahead))
      
      y.hat <- idx_series(
        as.matrix(rbind(y.t.t, y.hat.kfas)),
        start = positions[1])
      colnames(y.hat$data)<-c("y.hat","y.hat.upr","y.hat.lwr")
      
      i.level <- grep("level", colnames(att(model_output)))
      level.t.t <- idx_series(as.numeric(att(model_output)[, i.level]), start = positions[1])
      
      i.slope <- grep("slope", colnames(att(model_output)))
      slope.t.t <- idx_series(as.numeric(att(model_output)[, i.slope]), start = positions[1])
      
      if (!return.all) {
        cutoff <- tail(index, 1)
        keep <- positions > cutoff
        keep_positions <- positions[keep]
        if (length(keep_positions) > 0) {
          y.hat <- y.hat[keep_positions]
          level.t.t <- level.t.t[keep_positions]
          slope.t.t <- slope.t.t[keep_positions]
        }
      }
      
      out <- list(
        y.hat = y.hat,
        level.t.t = level.t.t,
        slope.t.t = slope.t.t,
        a.t.t = att(model_output),
        P.t.t = Ptt(model_output)
      )
      return(out)
    },
    get_growth_y = function(smoothed = FALSE, return.components = FALSE) {
      "Returns the growth rate of the incidence (\\eqn{y}) of the cumulated
      variable (\\eqn{Y}). Computed as
      \\deqn{g_t = \\exp\\{\\delta_t\\}+\\gamma_t.}
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{smoothed} Logical value indicating whether to use the
        smoothed estimates of \\eqn{\\delta} and \\eqn{\\gamma} to compute the
        growth rate (\\code{TRUE}), or the contemporaneous filtered estimates
        (\\code{FALSE}). Default is \\code{FALSE}.}
        \\item{\\code{return.components} Logical value indicating whether to
        return the estimates of \\eqn{\\delta} and \\eqn{\\gamma} as well as
        the estimates of the growth rate, or just the growth rate. Default is
        \\code{FALSE}.}
      }}
      \\subsection{Return Value}{\\code{idx_series} objects containing
      smoothed/filtered growth rates and components (\\eqn{\\delta} and
      \\eqn{\\gamma}), where applicable.}"
      kfs_out <- output
      idx <- index
      
      if (smoothed) {
        att <- alphahat(kfs_out)
      } else {
        att <- att(kfs_out)
      }
      
      filtered_slope <- idx_series(as.numeric(att[, "slope"]), start = idx[1])
      filtered.level <- idx_series(as.numeric(att[, "level"]), start = idx[1])
      g.t <- idx_series(exp(idx_values(filtered.level)), start = idx[1])
      gy.t <- idx_series(idx_values(g.t) + idx_values(filtered_slope), start = idx[1])
      if (return.components) {
        return(list(gy.t, g.t, filtered_slope))
      } else {
        return(gy.t)
      }
    },
    get_gy_ci = function(smoothed = FALSE, confidence.level = 0.68) {
      "Returns the growth rate of the incidence (\\eqn{y}) of the cumulated
      variable (\\eqn{Y}). Computed as
      \\deqn{g_t = \\exp\\{\\delta_t\\}+\\gamma_t.}
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{smoothed} Logical value indicating whether to use the
        smoothed estimates of \\eqn{\\delta} and \\eqn{\\gamma} to compute the
        growth rate (\\code{TRUE}), or the contemporaneous filtered estimates
        (\\code{FALSE}). Default is \\code{FALSE}.}
        \\item{\\code{confidence.level} Confidence level for the confidence
        interval.  Default is \\eqn{0.68}, which is one standard deviation for
        a normally distributed random variable.}
      }}
      \\subsection{Return Value}{\\code{idx_series} object containing smoothed/filtered
       growth rates and upper and lower bounds for the confidence intervals.}"
      
      kfs_out <- output
      idx <- index
      
      if (smoothed) {
        att <- alphahat(kfs_out)
        var <- get_V(kfs_out)
      } else {
        att <- att(kfs_out)
        var <- Ptt(kfs_out)
      }
      
      filtered_slope <- as.numeric(att[, "slope"])
      filtered.level <- as.numeric(att[, "level"])
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope
      
      idx.slope <- grep("slope", colnames(att(kfs_out)))
      ci <- qnorm((1 - confidence.level) / 2) *
        sqrt(as.numeric(var[idx.slope, idx.slope,])) %o% c(1, -1)
      ci_bounds <- gy.t + ci
      
      pred <- idx_series(cbind(fit = gy.t, lower=ci_bounds[,1], upper=ci_bounds[,2]), start = idx[1])
      
      return(pred)
    },
    print=function(){
      "Provides a quick glimpse of model states and standard errors."
      cat("Object of FilterResults Class\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
    },
    summary=function(){
      "Supplies details of the FilterResults object, such as estimated
      parameter values, start and end positions of estimation."
      H <- matrixKFS(output, "H")[, , 1]
      Q_gamma <- matrixKFS(output, "Q")[2, 2, 1]
      if (sea.period>1){  
        Q_seasonal <- matrixKFS(output, "Q")[3, 3, 1]
      }
      
      start.idx <- index[1]
      end.idx <- index[length(index)]
      
      cat("Summary of FilterResults Object\n")
      cat("Model Details:\n")
      cat("  - Estimation start position:", start.idx)
      cat("\n")
      cat("  - Estimation end position:", end.idx)
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
      if (ar1){
        ar1_comp<-matrixKFS(output,"T")["ar1","ar1",1]
        cat("  - AR(1) coefficient:", signif(ar1_comp,3))
        cat("\n")
      }
      cat("  - Variance parameter estimates\n")
      cat("Observation equation noise:",format(H, digits = 4))
      cat("\n")
      cat("State transition equation noise:",format(Q_gamma, digits = 4))
      cat("\n")
      cat("Signal-to-Noise Ratio (q):", format(Q_gamma / H, digits = 4))
      cat("\n")
      if (sea.period>1){
        cat("Seasonality noise:",format(Q_seasonal, digits = 4))
      }
    }, 
    mapes=function(n.ahead,Y){
      "Computes five metrics, including Mean Absolute Percentage Error (MAPE), 
      for forecasts against a holdout sample. For more details, please refer to 
    \\link{mapes}."
      if (xpred_logical){
        if (is.null(xpred.new)){
          stop("xpred.new cannot be NULL.")
        } 
      }
      p <- attr(modelKFS(output), 'p')
      if(p!=1) { stop('NotImplementedError') }
      
      estimation.end <- tail(index, 1)
      
      eval.window <- get_timeframe(Y, estimation.end, estimation.end + n.ahead)
      y.eval.diff <- idx_diff(eval.window, 1L)
      
      y.hat.diff.final <- .self$predict_level(
        n.ahead = n.ahead, confidence.level =0.68,
        sea.on = TRUE
      )
      
      # Extract the relevant, overlapping positions
      eval_pos <- idx_positions(y.eval.diff)
      eval_pos <- eval_pos[eval_pos > estimation.end]
      filtered_y_eval_diff <- y.eval.diff[eval_pos]
      forecast_mat <- idx_values(y.hat.diff.final)
      forecast_pos <- idx_positions(y.hat.diff.final)
      common_pos <- intersect(eval_pos, forecast_pos)
      
      d.eval <- data.frame(
        pos = common_pos,
        Actual = idx_values(filtered_y_eval_diff[common_pos]),
        Forecast = forecast_mat[match(common_pos, forecast_pos), 1],
        lwr = forecast_mat[match(common_pos, forecast_pos), 2],
        upr = forecast_mat[match(common_pos, forecast_pos), 3]
      )
      d.eval <- na.omit(d.eval)
      
      if (any(d.eval$Actual==0)){
        warning("Validation data contains zeros. MAPE is not a reliable measure.")
      }
      
      mape.sea <- mean(100*(abs(d.eval$Actual - d.eval$Forecast)/d.eval$Actual))
      smape<-mean(100*(abs(d.eval$Actual - d.eval$Forecast)/(d.eval$Actual+d.eval$Forecast)))
      mae<-abs(d.eval$Actual - d.eval$Forecast) %>% mean
      rmse<-sqrt(mean((d.eval$Actual - d.eval$Forecast)^2))
      coverage<-100*sum(d.eval$lwr<=d.eval$Actual & d.eval$upr>=d.eval$Actual)/n.ahead
      
      return(list(mape=mape.sea, smape=smape, mae=mae, rmse=rmse, coverage=coverage))
    }
  )
)