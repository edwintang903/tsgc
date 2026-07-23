# Refactored: analysis-only FilterResultsLI class operating on idx_series
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
#' @title Class for the estimated Leading Indicator Model
#'
#' @description A class that holds the information of an estimated Leading Indicator model.
#' Contains methods to extract smoothed/filtered estimates of the states, the
#' level of the incidence variable \eqn{y}, and forecasts of \eqn{y}. The output from the estimate method
#' of the SSModelLeadingIndicator class is of the class FilterResultsLI.
#'
#' @field data An \code{idx_series} object with cumulated variables: lagged leading
#' indicator and target variable (plus their increments and log-growth rates - see
#' \code{\link{add_daily_ldl}}).
#' @field output A \code{KFS} results object obtained after fitting a 
#' \code{SSModelLeadingIndicator} model.
#' @field n.lag Number of integer positions the leading indicator is lagged by, inherited
#' from the estimated \code{SSModelLeadingIndicator} model.
#' @field sea.period The period of seasonality, inherited from the estimated 
#' \code{SSModelLeadingIndicator} model. For a day-of-the-week
#'   effect with daily data, this would be 7. 
#' @field LeadIndCol The column in \code{data} that contains the leading 
#' indicator, inherited from the estimated \code{SSModelLeadingIndicator} model.
#' @field xpred_logical Vector of length 2 with logical values, indicating whether
#' there are exogenous predictors for leading series and target series respectively. 
#' @field xpred_lead.new An \code{idx_series} object containing the values of exogenous
#' variables for the leading indicator over the prediction time frame.
#' @field xpred_targ.new An \code{idx_series} object containing the values of exogenous
#' variables for the target variable over the prediction time frame.
#' @field start Integer position marking the start of the estimation period.
#' @field end Integer position marking the end of the estimation period.
#'
#' @references Harvey, A. (2021). TIME SERIES MODELLING OF EPIDEMICS:
#' LEADING INDICATORS, CONTROL GROUPS AND POLICY ASSESSMENT.
#' National Institute Economic Review, 257, 83-100.
#' doi:10.1017/nie.2021.21
#'
#' @importFrom magrittr %>%
#' @importFrom methods new setRefClass setOldClass
#' @importFrom stats predict qnorm
#' 
#' @examples
#' library(tsgc)
#' set.seed(1)
#' lead <- cumsum(rpois(150, 6)) + 1
#' targ <- cumsum(rpois(150, 8)) + 1
#' Y <- idx_series(cbind(lead, targ), start = 1)
#'
#' # Define and estimate the model
#' model <- SSModelLeadingIndicator(Y = Y, n.lag = 5, q = NULL, LeadIndCol = 1,
#'   sea.period = 0, start = 1, end = 100)
#' res <- estimate(model)
#'
#' # Print estimation results
#' res$print_estimation_results()
#'
#' # Forecast 7 positions ahead from the end of the estimation window
#' res$predict_level(n.ahead = 7, confidence.level = 0.68)
#'
#' # Forecast 7 positions ahead from the model and return filtered states
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
#' res$mapes(n.ahead = 7, Y = Y)
#'
#' @export
#'
FilterResultsLI <- setRefClass(
  "FilterResultsLI",
  fields = list(
    data = "idx_series",
    output = "KFS",
    n.lag="numeric",
    sea.period="numeric",
    LeadIndCol="numeric",
    xpred_lead.new="ANY",
    xpred_targ.new="ANY",
    xpred_logical="logical",
    start="ANY",
    end="ANY"),
  methods = list(
    initialize = function(data, output,n.lag,sea.period,LeadIndCol,
                          xpred_logical, start, end, 
                          xpred_lead.new=NULL, xpred_targ.new=NULL)
    {
      "Create an instance of the \\code{FilterResultsLI} class with fields defined
      earlier in the fields section."
      data <<- data
      output <<- output
      n.lag <<- n.lag
      sea.period <<- sea.period
      LeadIndCol<<-LeadIndCol
      start<<-start
      end<<-end
      xpred_lead.new<<-xpred_lead.new
      xpred_targ.new<<-xpred_targ.new
      xpred_logical<<-xpred_logical
    },
    predict_level = function(n.ahead=n.lag, 
                             confidence.level=0.68,
                             sea.on = TRUE){
      "Forecast the cumulated variable or the incidence of it. This function returns
      the forecast of the cumulated variable \\eqn{Y}, or the forecast of the incidence of the cumulated variable, \\eqn{y}. For
      example, in the case of an epidemic, \\eqn{y} might be daily new cases of
      the disease and
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of periods ahead you wish to forecast from
        the end of the estimation window. Default is \\code{n.lag}.}
        \\item{\\code{sea.on} Logical value indicating whether to return the prediction 
        of just the trend or prediction incorporating seasonality. Deafults to \\code{TRUE}.}
        \\item{\\code{confidence.level} The confidence level for the log growth
         rate that should be used to compute the forecast intervals of \\eqn{y}.}
       }
      }
      \\subsection{Return Value}{An \\code{idx_series} object containing the point
      forecasts and upper and lower bounds of the forecast interval.}"
      if (n.ahead==1){
        n.ahead=2
        unity=TRUE
      } else{
        unity=FALSE
      }
      
      if (!sea.on){
        # Create the forecasts
        # This gives the forecasts of delta
        forcout<-.self$predict_all(n.ahead, sea.on = FALSE, return.all = FALSE, 
                                   confidence.level=confidence.level)$y.hat.kfas
      } else {
        #Re-do with seasonal component
        forcout = .self$predict_all(n.ahead, sea.on = TRUE, return.all = FALSE)$y.hat.kfas
      }
      
      # Create empty matrix to put forecasts in
      forecasts <- matrix(NA_real_, ncol=3, nrow=n.ahead)
      colnames(forecasts) <- c('forc','lwr','upr')
      
      last_admit <- idx_values(data[end])[,"cTarg"]
      
      LDLtarg_fc <- forcout$LDLtarg
      
      # Compute forecasts as per (7) in Andrew's Time Series Models for Epidemics paper
      # Confidence intervals computed as per Harvey, Kattuman and Thamotheram 2021 NIESR paper
      forecasts[1, 1] = last_admit*exp(LDLtarg_fc[1,1])
      forecasts[2:n.ahead, 1] = last_admit*exp(LDLtarg_fc[2:n.ahead,1])*cumprod(1+exp(LDLtarg_fc[1:(n.ahead-1),1]))
      
      forecasts[1, 2] = last_admit*exp(LDLtarg_fc[1,2])
      forecasts[2:n.ahead, 2] = last_admit*exp(LDLtarg_fc[2:n.ahead,2])*cumprod(1+exp(LDLtarg_fc[1:(n.ahead-1),2]))
      
      forecasts[1, 3] = last_admit*exp(LDLtarg_fc[1,3])
      forecasts[2:n.ahead, 3] = last_admit*exp(LDLtarg_fc[2:n.ahead,3])*cumprod(1+exp(LDLtarg_fc[1:(n.ahead-1),3]))
      
      # Round forecasts to 2 decimal places
      forecasts <- round(forecasts, 2)
      
      fadmits <- idx_series(forecasts, start = end + 1L)
      
      if (unity){
        return(fadmits[end + 1L])
      } else{
        return(fadmits)
      }
    },
    print_estimation_results = function() {
      "Prints a table of estimated parameters in a format ready to paste into
      LaTeX."
      H1 <- output$model$H[1, 1, 1]
      H2 <- output$model$H[2, 2, 1]
      Q_gamma <- output$model$Q[2, 2, 1]
      has_seasonal <- sea.period > 1
      
      if (has_seasonal) {
        Q_seasonal <- output$model$Q[3, 3, 1]
        tbl <- data.frame(
          a = format(H1, digits = 3),
          b = format(H2, digits = 3),
          c = format(Q_gamma, digits = 3),
          d = format(Q_seasonal, digits = 3))
        header.names <- c('$\\sigma_\\varepsilon1^2$',
                          '$\\sigma_\\varepsilon2^2$',
                          '$\\sigma_{IRW}^2$',
                          '$\\sigma_{trend1}^2$')
      } else {
        tbl <- data.frame(
          a = format(H1, digits = 3),
          b = format(H2, digits = 3),
          c = format(Q_gamma, digits = 3))
        header.names <- c('$\\sigma_\\varepsilon1^2$',
                          '$\\sigma_\\varepsilon2^2$',
                          '$\\sigma_{IRW}^2$')
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
    predict_all = function(n.ahead, sea.on = TRUE, 
                           return.all = FALSE, 
                           confidence.level=0.68) {
      "Returns forecasts of the incidence variable \\eqn{y}, the state variables
       and the conditional covariance matrix for the states.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of forecasts you wish to create from
        the end of your sample period.}
        \\item{\\code{sea.on} Logical value indicating whether to return the prediction 
        of just the trend or prediction incorporating seasonality. Default is \\code{TRUE}.}
        \\item{\\code{return.all} Logical value indicating whether to return
        all filtered estimates and forecasts
        (\\code{TRUE}) or only the forecasts (\\code{FALSE}). Default is
        \\code{FALSE}.}
      }}
      \\subsection{Return Value}{\\code{idx_series} objects containing the forecast
      (and filtered, where applicable) level
      of \\eqn{y} (\\code{y.hat}), \\eqn{\\delta} (\\code{level.t.t}),
      \\eqn{\\gamma} (\\code{slope.t.t}), vector of states including the
      seasonals where applicable (\\code{a.t.t}) and covariance matrix of all
      states including seasonals where applicable (\\code{P.t.t}).}"
      new.model <- modelKFS(output)
      Qf = matrixKFS(output,"Q")[,,1]
      Hf = matrixKFS(output,"H")[,,1]
      oldn<-attr(new.model, 'n')
      
      # Provide observed leading indicator data
      na_vals<-matrix(NA, ncol = ncol(gety(new.model)), nrow = n.ahead)
      
      remaining_data <- idx_range(data)[2] - end
      future_rows<-min(n.ahead, n.lag, remaining_data)
      
      if (future_rows > 0) {
        true_leading <- get_timeframe(data, end + 1L, end + future_rows)
        na_vals[1:future_rows,1] = idx_values(true_leading)[,"LDLlead"]
      }
      
      #Supply the new data back to the new.model object
      new.model$y <- rbind(gety(new.model),na_vals) %>% as.ts()
      
      if (xpred_logical[1] || xpred_logical[2]){
        newZ<-array(new.model$Z[,,dim(new.model$Z)[3]], 
                    dim = c(dim(new.model$Z)[1], dim(new.model$Z)[2], n.ahead))
        if (xpred_logical[1]){
          if (is_idx_series(xpred_lead.new)){
            xpred_lead.new.subset<-as.matrix(idx_values(get_timeframe(
              xpred_lead.new, end + 1L, end + n.ahead)))
            d1<-ncol(xpred_lead.new.subset)
            newZ[1,1:d1,]<-t(xpred_lead.new.subset)
          } else {
            stop("xpred_lead.new not provided.")
          }
        }
        if (xpred_logical[2]){
          if (is_idx_series(xpred_targ.new)){
            xpred_targ.new.subset<-as.matrix(idx_values(get_timeframe(
              xpred_targ.new, end + 1L, end + n.ahead)))
            d2<-ncol(xpred_targ.new.subset)
            if (!xpred_logical[1]){d1=0}
            newZ[2,(d1+1):(d1+d2),]<-t(xpred_targ.new.subset)
          } else {
            stop("xpred_targ.new not provided.")
          }
        }
        new.model$Z <- abind::abind(new.model$Z,newZ,along = 3)
        attr(new.model, 'n') <- as.integer(oldn + n.ahead)
        model_output <- KFS(new.model)
        
        newdata <- if (sea.period<2 && !xpred_logical[1] && xpred_logical[2]){
          SSModel(na_vals ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+
                    SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1)+
                    SSMregression(~xpred_targ.new.subset, type="distinct", index=2),
                  H = Hf)
        } else if (sea.period<2 && xpred_logical[1] && !xpred_logical[2]){
          SSModel(na_vals ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+
                    SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1)+
                    SSMregression(~xpred_lead.new.subset, type="distinct", index=1),
                  H = Hf)
        } else if (sea.period<2 && xpred_logical[1] && xpred_logical[2]) {
          SSModel(na_vals ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+
                    SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1)+
                    SSMregression(~xpred_lead.new.subset, type="distinct", index=1)+
                    SSMregression(~xpred_targ.new.subset, type="distinct", index=2),
                  H = Hf)
        } else if (sea.period>=2 && !xpred_logical[1] && xpred_logical[2]) {
          SSModel(na_vals ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+
                    SSMseasonal(sea.period, Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                    SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1)+
                    SSMregression(~xpred_targ.new.subset, type="distinct", index=2),
                  H = Hf)
        } else if (sea.period>=2 && xpred_logical[1] && !xpred_logical[2]) {
          SSModel(na_vals ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+
                    SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                    SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1)+
                    SSMregression(~xpred_lead.new.subset, type="distinct", index=1),
                  H = Hf)
        } else {
          SSModel(na_vals ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+
                    SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                    SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1)+
                    SSMregression(~xpred_lead.new.subset, type="distinct", index=1)+
                    SSMregression(~xpred_targ.new.subset, type="distinct", index=2),
                  H = Hf)
        }
        
        if (sea.on == TRUE) {
          y.hat.kfas <- predict(
            modelKFS(output), interval = 'prediction',
            newdata = newdata, level = confidence.level, states = 'all')
        } else {
          y.hat.kfas <- predict(
            modelKFS(output), interval = 'prediction',
            newdata = newdata, level = confidence.level, states = 'level')
        }
        y.t.t<-matrix(nrow=2,ncol=oldn)
        
        for (j in 1:2){
          for (i in 1:oldn){
            y.t.t[j,i] <- output$att[i,] %*% drop(matrixKFS(output,"Z"))[j,,i]
          }
        }
        
      } else {
        attr(new.model, 'n') <- as.integer(oldn + n.ahead)
        model_output <- KFS(new.model)
        
        # Create forecast model object
        if (sea.period<2) {
          forcmodel = SSModel(na_vals ~ SSMtrend(degree = 2, 
                                                 Q = matrix(c(0,0,0,Qf[2,2]),2,2),
                                                 type = 'common')
                              +SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1),
                              H = matrixKFS(output,"H"))
        } else {
          forcmodel = SSModel(na_vals ~ SSMtrend(degree = 2, 
                                                 Q = matrix(c(0,0,0,Qf[2,2]),2,2),
                                                 type = 'common')
                              +SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), 
                                           sea.type='trigonometric', type='distinct')
                              +SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1),
                              H = matrixKFS(output,"H"))
        }
        
        if (sea.on == TRUE) {
          y.hat.kfas <- predict(
            output$model, interval = 'prediction',
            newdata = forcmodel, level = confidence.level, states = 'all')
        } else {
          y.hat.kfas <- predict(
            output$model, interval = 'prediction',
            newdata = forcmodel, level = confidence.level, states = 'level')
        }
        # Assumes time invariant Z.t
        y.t.t <- t(output$att %*% t(drop(matrixKFS(output,"Z"))))
      }
      
      n <- attr(output$model, "n")
      positions <- seq.int(start, length.out = (n + n.ahead))
      
      y.hat <- idx_series(
        c(y.t.t[2,], as.numeric(as.matrix(y.hat.kfas$LDLtarg[, 1]))),
        start = positions[1])
      
      i.level <- grep("level", colnames(att(model_output)))[1]
      level.t.t <- idx_series(as.numeric(att(model_output)[, i.level]), start = positions[1])
      i.slope <- grep("slope", colnames(att(model_output)))
      slope.t.t <- idx_series(as.numeric(att(model_output)[, i.slope]), start = positions[1])
      
      if (!return.all) {
        keep_positions <- positions[positions > end]
        if (length(keep_positions) > 0) {
          y.hat <- y.hat[keep_positions]
          level.t.t <- level.t.t[keep_positions]
          slope.t.t <- slope.t.t[keep_positions]
        }
      }
      
      out <- list(
        y.hat = y.hat,
        y.hat.kfas=y.hat.kfas,
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
      idx <- idx_positions(get_timeframe(data, start, end))
      
      if (smoothed) {
        att <- alphahat(output)
      } else {
        att <- att(output)
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
      
      idx <- idx_positions(get_timeframe(data, start, end))
      
      if (smoothed) {
        att <- alphahat(output)
        var <- get_V(output)
      } else {
        att <- att(output)
        var <- Ptt(output)
      }
      
      filtered_slope <- as.numeric(att[, "slope"])
      filtered.level <- as.numeric(att[, "level"])
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope
      
      idx.slope <- grep("slope", colnames(att(output)))
      ci <- qnorm((1 - confidence.level) / 2) *
        sqrt(as.numeric(var[idx.slope, idx.slope,])) %o% c(1, -1)
      ci_bounds <- gy.t + ci
      
      pred <- idx_series(cbind(fit=gy.t, lower=ci_bounds[,1], upper=ci_bounds[,2]), start = idx[1])
      
      return(pred)
    },
    print=function(){
      "Provides a quick glimpse of model states and standard errors."
      cat("Object of FilterResultsLI Class\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
    },
    summary=function(){
      "Supplies details of the FilterResultsLI object, such as estimated
      parameter values, start and end positions of estimation."
      H <- matrixKFS(output, "H")[, , 1]
      Q_gamma <- matrixKFS(output, "Q")[2, 2, 1]
      has_seasonal <- sea.period > 1
      cat("Summary of FilterResultsLI Object\n")
      cat("Model Details:\n")
      cat("  - Estimation start position:", start)
      cat("\n")
      cat("  - Estimation end position:", end)
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
      cat("  - Variance parameter estimates\n")
      cat("Observation equation noise:",format(H, digits = 4))
      cat("\n")
      cat("State transition equation noise:",format(Q_gamma, digits = 4))
      if (has_seasonal) {
        Q_seasonal <- matrixKFS(output, "Q")[3, 3, 1]
        cat("\n")
        cat("Seasonality noise:",format(Q_seasonal, digits = 4))
      }
    },
    mapes=function(n.ahead,Y){
      "Computes five metrics, including Mean Absolute Percentage Error (MAPE), 
      for forecasts against a holdout sample. For more details, please refer to 
    \\link{mapes}."
      sea<-.self$predict_level(n.ahead=n.ahead, sea.on=TRUE)
      
      eval_window <- get_timeframe(Y, end)
      data_validation <- add_daily_ldl(eval_window, LeadIndCol = LeadIndCol)
      newTarg_validation <- data_validation$newTarg
      
      common_pos <- intersect(idx_positions(newTarg_validation), idx_positions(sea))
      common_pos <- head(common_pos, n.ahead)
      if (length(common_pos) == 0) {
        stop("No overlapping positions between the holdout sample and the forecast horizon.")
      }
      
      sea_mat <- idx_values(sea)
      d.eval <- data.frame(
        Actual = idx_values(newTarg_validation[common_pos]),
        Forecast = sea_mat[match(common_pos, idx_positions(sea)), 1],
        lwr = sea_mat[match(common_pos, idx_positions(sea)), 2],
        upr = sea_mat[match(common_pos, idx_positions(sea)), 3]
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