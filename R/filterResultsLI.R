setOldClass("KFS")
#'
#' @title Class for the estimated Leading Indicator Model
#'
#' @description A class that holds the information of an estimated Leading Indicator model.
#' Contains methods to extract smoothed/filtered estimates of the states, the
#' level of the incidence variable \eqn{y}, and forecasts of \eqn{y}. The output from the estimate method
#' of the SSModelLeadingIndicator class is of the class FilterResultsLI.
#'
#' @field data_xts A xts object with cumulated variables: lagged leading indicator and 
#' target variable.
#' @field output A \code{KFS} results object obtained after fitting a 
#' \code{SSModelLeadingIndicator} model.
#' @field n.lag Number of days to lag the leading indicator, inherited from the 
#' estimated \code{SSModelLeadingIndicator} model.
#' @field sea.period The period of seasonality, inherited from the estimated 
#' \code{SSModelLeadingIndicator} model. For a day-of-the-week
#'   effect with daily data, this would be 7. 
#' @field LeadIndCol The column in \code{data_xts} that contains the leading 
#' indicator, inherited from the estimated \code{SSModelLeadingIndicator} model.
#' @field xpred_logical Vector of length 2 with logical values, indicating whether
#' there are exogenous predictors for leading series and target series. 
#' @field xpred_lead.new An xts object containing the values of exogenous variables for 
#' the leading indicator over the prediction time frame.
#' @field xpred_targ.new An xts object containing the values of exogenous variables for 
#' the target variable over the prediction time frame.
#' @field resolution A character object showing the time resolution of the data 
#' in \code{data_xts}. Options are "daily", "monthly, "quarterly" and "yearly".
#' Automatically estimated when \code{data_xts} is provided.
#'
#'@references Harvey, A. (2021). TIME SERIES MODELLING OF EPIDEMICS: 
#'LEADING INDICATORS, CONTROL GROUPS AND POLICY ASSESSMENT. 
#'National Institute Economic Review, 257, 83–100. 
#'doi:10.1017/nie.2021.21
#'
#' @importFrom xts periodicity last lag.xts
#' @importFrom magrittr %>%
#' @importFrom methods new
#' 
#' @examples
#' library(tsgc)
#' 
#' #Define estimation environment
#' Y = england[,1:2] #select only cumulative cases and cumulative hospitalizations
#' estimation.date.start = as.Date("2021-04-30")
#' estimation.date.end = as.Date("2021-07-24")
#' plt.length=30
#' q = NULL
#' confidence.level = 0.68
#' plt.length = 14
#' n.lag = 4
#' 
#' # Define and Estimate the model
#' out<-SSModelLeadingIndicator(Y=Y,n.lag = n.lag,q=q,LeadIndCol=1, sea.period=7,
#' start.date=estimation.date.start, end.date=estimation.date.end) 
#' res<-estimate(out) 
#' 
#' # Print estimation results
#' res$print_estimation_results()
#' 
#' # Forecast 7 days ahead from the end of the estimation window
#' res$predict_level(n.ahead = 7, confidence.level = 0.68)
#'   
#' # Forecast 7 days ahead from the model and return filtered states
#' res$predict_all(n.ahead = 7, return.all = TRUE)
#' 
#' # Return the filtered growth rate and its components
#' res$get_growth_y(return.components = TRUE)
#' 
#' # Return smoothed growth rate of incidence variable and its confidence
#' # interval
#' res$get_gy_ci(smoothed = TRUE, confidence.level = 0.68)
#'
#' # Plot forecast of new cases 7 days ahead
#' res$plot_log_forecast(Y=Y, n.ahead=7,
#' plt.start.date=estimation.date.end-plt.length)
#' 
#' # Plot forecast of new cases 7 days ahead
#' res$plot_forecast(n.ahead=7,
#' plt.start.date = estimation.date.end-plt.length,
#' series.name="hospitalizations")
#' 
#' # Plot forecasts and outcomes over evaluation period
#' res$plot_holdout(Y=Y,n.ahead=7, series.name="hospitalizations")
#' 
#' # Plot filtered gy, g and gamma
#' res$plot_gy_components(plt.start.date = estimation.date.end-plt.length)
#' 
#' # Plot filtered gy, g and gamma
#' res$plot_gy_ci(plt.start.date = estimation.date.end-plt.length)
#' 
#' #Return MAPE of forecast
#' res$mapes(n.ahead=7,Y)
#'
#' @export
#'
FilterResultsLI <- setRefClass(
  "FilterResultsLI",
  fields = list(
    data_xts = "ANY",
    output = "KFS",
    n.lag="numeric",
    sea.period="numeric",
    LeadIndCol="numeric",
    xpred_lead.new="ANY",
    xpred_targ.new="ANY",
    xpred_logical="logical",
    resolution="character",
    start.date="ANY",
    end.date="ANY"),
  methods = list(
    initialize = function(data_xts, output,n.lag,sea.period,LeadIndCol,
                          xpred_logical, start.date, end.date, 
                          xpred_lead.new=NULL, xpred_targ.new=NULL, resolution="daily")
    {
      "Create an instance of the \\code{FilterResultsLI} class with fields defined
      earlier in the fields section."
      data_xts <<- data_xts
      output <<- output
      n.lag <<- n.lag
      sea.period <<- sea.period
      LeadIndCol<<-LeadIndCol
      start.date<<-start.date
      end.date<<-end.date
      xpred_lead.new<<-xpred_lead.new
      xpred_targ.new<<-xpred_targ.new
      xpred_logical<<-xpred_logical
      resolution<<-get_time_resolution(index(data_xts))
    },
    predict_level = function(n.ahead=n.lag, 
                             confidence.level=0.68,
                             sea.on = FALSE){
      "Forecast the cumulated variable or the incidence of it. This function returns
      the forecast of the cumulated variable \\eqn{Y}, or the forecast of the incidence of the cumulated variable, \\eqn{y}. For
      example, in the case of an epidemic, \\eqn{y} might be daily new cases of
      the disease and
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of periods ahead you wish to forecast from
        the end of the estimation window. Default is \\code{n.lag}.}
        \\item{\\code{sea.on} Logical value indicating whether to return the prediction 
        of just the trend or prediction incorporating seasonality.}
        \\item{\\code{confidence.level} The confidence level for the log growth
         rate that should be used to compute the forecast intervals of \\eqn{y}.}
       }
      }
      \\subsection{Return Value}{A list object containing n.lag and 2 \\code{xts}
      objects: the point forecasts and upper and lower bounds of the forecast interval
      for trend and forecast with seasonal component.}"
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
      
      # Create empty dataframe to put forecasts in
      forecasts <- matrix(NA,ncol=2,nrow=n.ahead) %>%
        as.data.frame()
      colnames(forecasts) = c('Admissions','Cases')
      
      last_cases<-as.vector(get_timeframe(data_xts, end.date)[1,]$cLead)
      last_admit<-as.vector(get_timeframe(data_xts, end.date)[1,]$cTarg)
      
      # Compute forecasts as per (7) in Andrew's Time Series Models for Epidemics paper
      # Confidence intervals computed as per Harvey, Kattuman and Thamotheram 2021 NIESR paper
      #forecasts$Cases[1] = last_cases*exp(forcout$LDLcases[1,1])
      #forecasts$Cases[2:n.ahead] = last_cases*exp(forcout$LDLcases[2:n.ahead,1])*cumprod(1+exp(forcout$LDLcases[1:(n.ahead-1),1]))
      #forecasts$Cases.lwr[1] = last_cases*exp(forcout$LDLcases[1,2])
      #forecasts$Cases.lwr[2:n.ahead] = last_cases*exp(forcout$LDLcases[2:n.ahead,2])*cumprod(1+exp(forcout$LDLcases[1:(n.ahead-1),2]))
      #forecasts$Cases.upr[1] = last_cases*exp(forcout$LDLcases[1,3])
      #forecasts$Cases.upr[2:n.ahead] = last_cases*exp(forcout$LDLcases[2:n.ahead,3])*cumprod(1+exp(forcout$LDLcases[1:(n.ahead-1),3]))
      
      forecasts$Admissions[1] = last_admit*exp(forcout$LDLtarg[1,1])
      forecasts$Admissions[2:n.ahead] = last_admit*exp(forcout$LDLtarg[2:n.ahead,1])*cumprod(1+exp(forcout$LDLtarg[1:(n.ahead-1),1]))
      
      forecasts$Admissions.lwr[1] = last_admit*exp(forcout$LDLtarg[1,2])
      forecasts$Admissions.lwr[2:n.ahead] = last_admit*exp(forcout$LDLtarg[2:n.ahead,2])*cumprod(1+exp(forcout$LDLtarg[1:(n.ahead-1),2]))
      
      forecasts$Admissions.upr[1] = last_admit*exp(forcout$LDLtarg[1,3])
      forecasts$Admissions.upr[2:n.ahead] = last_admit*exp(forcout$LDLtarg[2:n.ahead,3])*cumprod(1+exp(forcout$LDLtarg[1:(n.ahead-1),3]))
      
      # Round forecasts to nearest whole number
      forecasts = round(forecasts,2)
      
      # Put forecasts into a separate dataframe for admissions and cases
      admissions_forecasts = cbind(forecasts$Admissions,forecasts$Admissions.lwr,forecasts$Admissions.upr)
      colnames(admissions_forecasts) = c('forc','lwr','upr')
      
      #cases_forecasts = cbind(forecasts$Cases,forecasts$Cases.lwr,forecasts$Cases.upr)
      #colnames(cases_forecasts) = c('forc','lwr','upr')
      
      #Save forecast dates
      startforc = if (resolution=="daily"||resolution=="yearly"){
        end.date+1
      } else if (resolution=="quarterly"){
        end.date+0.25
      } else if (resolution=="monthly"){
        end.date+1/12
      }
      
      finds = seq_dates(startforc, length.out=n.ahead, resolution=resolution)
      
      fadmits = xts(admissions_forecasts[1:n.ahead,],finds)
      
      if (unity){
        return(fadmits[1,])
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
    predict_all = function(n.ahead, sea.on = FALSE, return.all = FALSE, 
                           confidence.level=0.68) {
      "Returns forecasts of the incidence variable \\eqn{y}, the state variables
       and the conditional covariance matrix
      for the states.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of forecasts you wish to create from
        the end of your sample period.}
        \\item{\\code{sea.on} Logical value indicating whether seasonal
        components should be included in the
        state-space model or not. Default is \\code{FALSE}.}
        \\item{\\code{return.all} Logical value indicating whether to return
        all filtered estimates and forecasts
        (\\code{TRUE}) or only the forecasts (\\code{FALSE}). Default is
        \\code{FALSE}.}
      }}
      \\subsection{Return Value}{\\code{xts} object containing the forecast
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
      
      new_index<-which.max(index(data_xts)==end.date)
      remaining_data<-dim(data_xts)[1]-new_index
      future_rows<-min(n.ahead, n.lag, remaining_data)
      last_index<-new_index+future_rows
      
      true_leading<-data_xts[(new_index+1):last_index]
      
      na_vals[1:future_rows,1] = as.vector(true_leading$LDLlead)
      
      #Supply the new data back to the new.model object
      new.model$y <- rbind(gety(new.model),na_vals) %>% as.ts()
      
      if (xpred_logical[1] || xpred_logical[2]){
        newZ<-array(new.model$Z[,,dim(new.model$Z)[3]], 
                    dim = c(dim(new.model$Z)[1], dim(new.model$Z)[2], n.ahead))
        if (xpred_logical[1]){
          if (is.xts(xpred_lead.new)){
            new_index<-which.max(index(data_xts)==end.date)
            xpred_lead.new.subset<-lag(xpred_lead.new,n.lag)[(new_index+1):(new_index+n.ahead)]
            d1<-dim(xpred_lead.new.subset)[2]
            newZ[1,1:d1,]<-t(xpred_lead.new.subset)
          } else {
            stop("xpred_lead.new not provided.")
          }
        }
        if (xpred_logical[2]){
          if (is.xts(xpred_targ.new)){
            new_index<-which.max(index(data_xts)==end.date)
            xpred_targ.new.subset<-xpred_targ.new[(new_index+1):(new_index+n.ahead)]
            d2<-dim(xpred_targ.new.subset)[2]
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
        if (sea.period==0) {
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
      
      dates <- seq_dates(start.date, length.out = (n + n.ahead), 
                         resolution=resolution)

      y.hat <- xts::xts(
        c(y.t.t[2,], y.hat.kfas$LDLtarg[, 1] %>% as.matrix()),
        order.by = dates)

      i.level <- grep("level", colnames(att(model_output)))[1]
      level.t.t <- xts::xts(att(model_output)[, i.level], order.by = dates) %>%
        as.xts()
      i.slope <- grep("slope", colnames(att(model_output)))
      slope.t.t <- xts::xts(att(model_output)[, i.slope], order.by = dates) %>%
        as.xts()
      
      newdate<-seq_dates(end.date,resolution=resolution,length.out=2)[2]
      if (!return.all) {
        y.hat <- get_timeframe(y.hat, newdate)
        level.t.t <- get_timeframe(level.t.t, newdate)
        slope.t.t <-  get_timeframe(slope.t.t, newdate)
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
      \\subsection{Return Value}{\\code{xts} object containing
      smoothed/filtered growth rates and components (\\eqn{\\delta} and
      \\eqn{\\gamma}), where applicable.}"
      idx <- index(get_timeframe(data_xts, start.date, end.date))

      if (smoothed) {
        att <- alphahat(output)
      } else {
        att <- att(output)
      }

      filtered_slope <- xts(att[, "slope"], order.by = idx)
      filtered.level <- xts(att[, "level"], order.by = idx)
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope
      names(gy.t) <- if (smoothed) { "smoothed gy.t" } else { "filtered gy.t" }
      names(g.t) <- if (smoothed) { "smoothed g.t" } else { "filtered g.t" }
      names(filtered_slope) <- if (smoothed) { "smoothed gamma.t" } else {
        "filtered gamma.t" }
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
      \\subsection{Return Value}{\\code{xts} object containing smoothed/filtered
       growth rates and upper and lower bounds for the confidence intervals.}"
      
      idx <- index(get_timeframe(data_xts, start.date, end.date))
      
      if (smoothed) {
        att <- alphahat(output)
      } else {
        att <- att(output)
      }
      
      filtered_slope <- xts(att[, "slope"], order.by = idx)
      filtered.level <- xts(att[, "level"], order.by = idx)
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope
      
      idx.slope <- grep("slope", colnames(att(output)))
      ci <- qnorm((1 - confidence.level) / 2) *
        sqrt(output$Ptt[idx.slope, idx.slope,]) %o% c(1, -1)
      ci_bounds <- as.vector(gy.t) + ci
      
      pred <- xts(cbind(gy.t, ci_bounds), order.by = idx)
      colnames(pred) <- c("fit","lower","upper")
      
      return(pred)
    },
    print=function(){
      "Provides a quick glimpse of model states and standard errors."
      cat("Object of FilterResultsLI Class\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
    },
    summary=function(){
      "Supplies details of the filterResults object, such as estimated 
      parameter values, start and end dates of estimation."
      H <- matrixKFS(output, "H")[, , 1]
      Q_gamma <- matrixKFS(output, "Q")[2, 2, 1]
      Q_seasonal <- matrixKFS(output, "Q")[3, 3, 1]
      cat("Summary of FilterResultsLI Object\n")
      cat("Model Details:\n")
      if (resolution=="daily"){
        cat("  - Estimation start date:", format(as.Date(start.date, origin = "1970-01-01"))) 
        cat("\n")
        cat("  - Estimation end date:", format(as.Date(end.date, origin = "1970-01-01")))
      } else if (resolution=="quarterly"){
        cat("  - Estimation start date:", format(as.yearqtr(start.date))) 
        cat("\n")
        cat("  - Estimation end date:", format(as.yearqtr(end.date)))
      } else if (resolution=="monthly"  || resolution=="yearly"){
        cat("  - Estimation start date:", format(as.yearmon(start.date))) 
        cat("\n")
        cat("  - Estimation end date:", format(as.yearmon(end.date)))
      } 
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
      cat("  - Variance parameter estimates\n")
      cat("Observation equation noise:",format(H, digits = 4))
      cat("\n")
      cat("State transition equation noise:",format(Q_gamma, digits = 4))
      cat("\n")
      cat("Seasonality noise:",format(Q_seasonal, digits = 4))
    },
    plot_forecast = function(n.ahead=7, confidence.level = 0.68, 
                              title=NULL, plt.start.date=NULL, 
                              series.name="target variable")
    {
      "Generates a forecast plot for the difference in the cumulative target 
      variable, showing actual values, forecasts including seasonal components,
      and prediction intervals around the forecasts. For more details, see 
      \\link{plot_forecast}."
        
        if (is.null(plt.start.date)){plt.start.date <- start.date}
        # add forecasts to plotting dataframe
        # fadmits<-.self$predict_level(n.ahead=n.ahead, 
        #                              confidence.level=confidence.level, 
        #                              sea.on=FALSE)
        sea<-.self$predict_level(n.ahead=n.ahead, 
                                 confidence.level=confidence.level, 
                                 sea.on=TRUE)
        
        # Create smoothed admissions
        # lcadmit = lag(data_xts$cAdmit) %>% na.omit()
        #smldlh = predict(output$model,states='trend')$LDLhosp %>% exp %>% as.vector
        #smadmit = smldlh*lcadmit
        #smAdmit = smadmit %>% xts(index(data_xts[(n.lag+1):(length(lcadmit)),])+1)
        
        #Plot forecast graph
        df_plot<-get_timeframe(data_xts, start.date, end.date)$newTarg
        #df_plot$Smooth<-smAdmit
        df_plot$Forecast<-sea[,1]
        #df_plot$ForecastTrend<-fadmits$forc
        df_plot<-get_timeframe(df_plot,plt.start.date)
        dates<-if (resolution=='quarterly'){
          qtr2date(as.yearqtr(index(df_plot)))
        } else if (resolution=='monthly' || resolution=="yearly"){
          qtr2date(as.yearmon(index(df_plot)))
        } else {
          as.Date(index(df_plot))
        }
        
        df_plot<-as.data.frame(df_plot)
        df_plot$Date <- dates
        
        ci_plot<-sea
        dates2<-if (resolution=='quarterly'){
          qtr2date(as.yearqtr(index(ci_plot)))
        } else if (resolution=='monthly' || resolution=="yearly"){
          qtr2date(as.yearmon(index(ci_plot)))
        } else {
          as.Date(index(ci_plot))
        }
        ci_plot<-as.data.frame(ci_plot)
        ci_plot$Date <- dates2
        
        p2<-ggplot2::ggplot(data = df_plot, aes(x = Date)) +
          ggplot2::geom_line(aes(y = newTarg, color = "Data"), lwd = 0.85) +
          #ggplot2::geom_line(aes(y = Smooth, color = "Smoothed\ndata"),lwd=0.85)+
          ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
          # ggplot2::geom_line(
          #   aes(y = ForecastTrend, color = "Forecast\nTrend"), lwd = 0.85
          # ) +
          ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
          ggplot2::geom_ribbon(data = ci_plot, aes(x = Date, ymin = lwr, ymax = upr),
                               linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
          labs(x = "Date", y = paste("New",series.name), title = title) +
          theme_economist_white(gray_bg = FALSE, base_size = 12) +
          theme(legend.title = element_blank()) +
          theme(
            text = element_text(size = rel(1.1)),
            axis.text = element_text(size = rel(1)),
            axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
            axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
            plot.title = element_text(margin=margin(b=5)),
            plot.caption = element_text(size = rel(1))
          ) +
          ggplot2::scale_linetype_manual(
            values = c("solid", "solid")) +
          ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
          ggplot2::scale_size_manual(values = c(1, 1, 1))
        return(p2)
      },
  plot_log_forecast = function(Y, n.ahead = 14,
                               plt.start.date=NULL, title="", caption = ""){
    "Plots actual and filtered values of the log cumulative growth rate 
      (\\eqn{\\ln(g_t)}) of the target variable in the estimation sample and 
      the forecast and realised log cumulative growth rate of the target variable
      out of the estimation sample. For more details, see \\link{plot_log_forecast}."
    
    forcout_sea<-.self$predict_all(n.ahead, sea.on = TRUE, return.all = FALSE)$y.hat
    old<-get_timeframe(data_xts, start.date, end.date)[,"LDLtarg"]
    
    eng_full<-add_daily_ldl(Y, LeadIndCol = LeadIndCol)
    eng_full<-eng_full[index(eng_full)>end.date,"LDLtarg"]
    
    actual<-eng_full[1:min(n.ahead, dim(eng_full)[1])]
    
    # Show filtered level only when xpred_logical is both FALSE
    if (!any(xpred_logical)){
      forcout<-.self$predict_all(n.ahead, sea.on = FALSE, return.all = FALSE)$y.hat
      smldlh = predict(output$model,states='trend')$LDLtarg
      
      dates_filtered <- seq_dates(from = start.date, length.out = length(smldlh), resolution=resolution)
      
      filtered=as.xts(as.vector(smldlh), order.by=dates_filtered)
      d.plot<-cbind(old,rbind(filtered,forcout),forcout_sea,actual)
      colnames(d.plot)<-c('EstimationSample', 'FilteredLevel', 'Forecast', 'RealisedData')
      
      #Graphical settings
      linetype_values<-c("solid", "solid", "solid", "dashed")
      color_values<-c(1, 2, 3, 'grey')
    } else {
      d.plot<-cbind(old,forcout_sea,actual)
      colnames(d.plot)<-c('EstimationSample', 'Forecast', 'RealisedData')
      
      #Graphical settings
      linetype_values<-c("solid", "solid", "dashed")
      color_values<-c(1, 3, 'grey')
    }
    
    if (!is.null(plt.start.date)){
      d.plot<-d.plot[index(d.plot)>=plt.start.date]
      df_plot <- as.data.frame(d.plot)
    } else{
      df_plot <- as.data.frame(d.plot)
    }
    df_plot$Date <- if (resolution=='quarterly'){
      qtr2date(as.yearqtr(rownames(df_plot)))
    } else if (resolution=='monthly' || resolution=="yearly"){
      qtr2date(as.yearmon(rownames(df_plot)))
    } else {
      as.Date(rownames(df_plot))
    }
    
    p1<-ggplot2::ggplot(data = df_plot, aes(x = Date)) +
      ggplot2::geom_line(aes(
        y = EstimationSample, color = "Estimation\nSample"), lwd = 0.85)
    
    if (!any(xpred_logical)){
      p1<-p1 +ggplot2::geom_line(aes(y = FilteredLevel, color = "Filtered\nLevel"),
                           lwd = 0.85)}
    
    p1<-p1+
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData"),
                         lwd = 0.85) +
      ggplot2::scale_color_manual(values = color_values) +
      scale_linetype_manual(
        values = linetype_values) +
      scale_x_date(labels = scales::date_format("%d %b %y")) +
      labs(x = "Date", y = "Log Growth Rate", caption = caption,
           title = title
      ) +
      theme_economist_white(gray_bg = FALSE) +
      scale_fill_economist() +
      theme(legend.title = element_blank()) +
      theme(
        text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
        axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
        plot.title = element_text(margin=margin(b=5)),
        plot.caption = element_text(size = rel(1)))
    return(p1)
  }, 
  plot_gy_components = function(plt.start.date = NULL,
                                 smoothed = FALSE, title = NULL){
    "Plots the growth rates and slope of the log cumulative growth rate 
      against the dates in estimation sample. 
      For more details, please see \\link{plot_gy_components}."
    Date <- Value <- Variable <- NULL
    # Determine plot start date
    if(is.null(plt.start.date)) {
        plt.start.date <- start.date
    }
    
    # Get gy.t, g.t and gamma
    gy.components <- .self$get_growth_y(return.components = TRUE, smoothed =
                                        smoothed)
    gy.t <- gy.components[[1]]
    g.t <- gy.components[[2]]
    gamma.t <- gy.components[[3]]
    
    d <- cbind(gy.t,g.t,gamma.t)
    names(d) <- c('gy.t','g.t','gamma.t')
    
    df_plot <- as.data.frame(d)
    df_plot$Date <- if (resolution=='quarterly'){
      qtr2date(as.yearqtr(rownames(df_plot)))
    } else if (resolution=='monthly' || resolution=="yearly"){
      qtr2date(as.yearmon(rownames(df_plot)))
    } else {
      as.Date(rownames(df_plot))
    }
    
    df_long <- df_plot %>%
      dplyr::filter(Date >= plt.start.date) %>%
      pivot_longer(cols = c(gy.t, g.t, gamma.t), names_to = "Variable",
                   values_to = "Value")
    
    p1 <- ggplot(df_long, aes(x = Date, y = Value, color = Variable)) +
      geom_line(lwd=0.85) +
      ggplot2::facet_wrap(~ factor(
        Variable, c("gy.t", "g.t", "gamma.t")), ncol = 1, scales = "free_y") +
      labs(title = title, y=ggplot2::element_blank()) +
      scale_color_manual(values = c("#AA2045","darkgrey","black")) +
      scale_x_date(labels = scales::date_format("%d %b %y")) +
      scale_y_continuous(breaks = waiver(), n.breaks = 4) +
      theme_economist_white(gray_bg = FALSE, base_size = 14) +
      theme(text = element_text(size= rel(1), margin=ggplot2::margin(b=5)),
            axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
            legend.position = "none")
    
    return(p1)
  },
  plot_gy_ci =function(plt.start.date = NULL, smoothed = FALSE,
                           title = NULL, series.name = NULL, pad.right = NULL){
    "Plots the growth rates and the slope of the log cumulative growth rate of 
    the target variable against the dates in estimation sample. 
      For more details, please see \\link{plot_gy_ci}."
    Date <- fit <- upper <- lower <- NULL
    
    # Determine plot start date
    if(is.null(plt.start.date)) {
        plt.start.date <- start.date
    }
    
    # Get confidence intervals to plot
    gy.ci<- .self$get_gy_ci(smoothed = smoothed)
    
    y.lab <- if(is.null(series.name)) { c("Growth rate") } else {
      paste("Growth rate of"," ",series.name,sep="")
    }
    
    df_plot <- as.data.frame(gy.ci)
    df_plot$Date <- if (resolution=='quarterly'){
      qtr2date(as.yearqtr(rownames(df_plot)))
    } else if (resolution=='monthly' || resolution=="yearly"){
      qtr2date(as.yearmon(rownames(df_plot)))
    } else {
      as.Date(rownames(df_plot))
    }
    
    p1 <- ggplot2::ggplot(df_plot[df_plot$Date>=plt.start.date,], aes(x=Date)) +
      ggplot2::geom_line(aes(y = fit), lwd = 0.85) +
      ggplot2::geom_hline(yintercept=0, linetype="solid",
                          color = "green", linewidth=1)+
      ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper),
                           linetype = 0, linewidth = 0, fill = "#AA2045",
                           alpha = 0.3) +
      ggplot2::scale_color_manual(values = c("black")) +
      geom_hline(
        aes(yintercept = 0.0), linetype = "solid", color = "green", lwd = 1.
      ) +
      labs(title=title, x="Date", y=y.lab) +
      theme_economist_white(gray_bg = FALSE, base_size = 14) +
      theme(
        legend.title = element_blank(),
        text = element_text(size = rel(1.)),
        axis.text = element_text(size = rel(1.)),
        axis.title.y = element_text(
          size = rel(1.),margin = ggplot2::margin(r=10)),
        axis.title.x = element_text(
          size = rel(1.),margin = ggplot2::margin(t=10)),
        plot.caption = element_text(size = rel(1))
      ) +
      theme(panel.grid.major.x = ggplot2::element_line(
        color = "gray50", linewidth = 0.5)) +
      scale_linetype_manual(
        values = c("solid")) +
      scale_x_date(labels = scales::date_format("%d %b %y"))
    
    if (!is.null(pad.right)) {
      plt.end.date <- tail(index(gy.ci),1)
      p1 <- p1 +
        ggplot2::scale_x_date(
          limits = c(as.Date(plt.start.date), plt.end.date + pad.right))
    }
    
    return(p1)
  },
  plot_holdout=function(Y,n.ahead=14, confidence.level = 0.68,
                        series.name = "target variable",
                        title= NULL, caption = NULL){
    "Plots the forecast of the target variable over a holdout sample. 
    For more details, please refer to \\link{plot_holdout}."
    # fadmits<-.self$predict_level(n.ahead=n.ahead, 
    #                              confidence.level=confidence.level, 
    #                              sea.on=FALSE)
    sea<-.self$predict_level(n.ahead=n.ahead, 
                             confidence.level=confidence.level, 
                             sea.on=TRUE)
    
    newdate<-seq_dates(end.date,resolution=resolution,length.out=2)[2]
    future_data<-get_timeframe(add_daily_ldl(Y,LeadIndCol = .self$LeadIndCol), newdate)
    if (n.ahead>dim(future_data)[1]){
      stop("The number of entries in the holdout sample is shorter than n.ahead.
           Please choose a smaller n.ahead, shorten the estimation period or 
           provide more holdout data.")
    } else {
      data_validation<-future_data[1:n.ahead, c("cTarg", "newTarg")]
    }
    
    newTarg_validation<-data_validation[,c("newTarg")]
    compare<-cbind(newTarg_validation, sea[,1])
    names(compare)<-c("Actual", "Forecast")
    
    # mape.trend <- 100*(abs(compare$Actual - compare$ForecastTrend)/
    #                      compare$Actual) %>% mean %>% round(4)
    mape.sea <- 100*(abs(compare$Actual - compare$Forecast)/compare$Actual) %>%
      mean %>% round(4)
    
    ci<-sea[,-1]
    colnames(ci) <- c('lower', 'upper')
    
    df_plot <- as.data.frame(compare)
    df_plot$Date <- if (resolution=='quarterly'){
      qtr2date(as.yearqtr(rownames(df_plot)))
    } else if (resolution=='monthly' || resolution=="yearly"){
      qtr2date(as.yearmon(rownames(df_plot)))
    } else {
      as.Date(rownames(df_plot))
    }
    
    ci_plot <- as.data.frame(ci)
    ci_plot$Date <- if (resolution=='quarterly'){
      qtr2date(as.yearqtr(rownames(ci_plot)))
    } else if (resolution=='monthly' || resolution=="yearly"){
      qtr2date(as.yearmon(rownames(ci_plot)))
    } else {
      as.Date(rownames(ci_plot))
    }
    
    
    p1<-ggplot2::ggplot(data = df_plot, aes(x = Date)) +
      ggplot2::geom_line(aes(y = Actual, color = "Actual"),lwd = 0.85) +
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast"),lwd = 0.85) +
      # ggplot2::geom_line(
      #   aes(y = ForecastTrend, color = "Forecast\nTrend"),lwd = 0.85) +
      ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
      ggplot2::geom_ribbon(data = ci_plot, aes(x = Date, ymin = lower, ymax = upper),linetype = 0, linewidth = 0, fill = "#AA2045",
                           alpha = 0.1) +
      labs(x = "Date", y = paste("New",series.name), title = title,
           subtitle = paste("MAPE: ",mape.sea,"%.",sep="")) +
      theme_economist_white(gray_bg = FALSE, base_size = 14) +
      theme(legend.title = element_blank()) +
      theme(
        text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1), margin = margin(r=10)),
        axis.title.x = element_text(size = rel(1), margin = margin(t=10)),
        plot.title = element_text(margin=margin(b=5)),
        plot.subtitle = element_text(
          size = rel(1), hjust=0,  margin = margin(t=3))
      ) +
      scale_linetype_manual(
        values = c("solid", "solid")) +
      scale_x_date(labels = scales::date_format("%d %b %y")) +
      scale_size_manual(values = c(1, 1.5, 1))
    
    return(p1)
  },
  mapes=function(n.ahead,Y){
    "Compute Mean Absolute Percentage Error (MAPE) for trend and seasonal 
    forecasts against a holdout sample. For more details, please refer to 
    \\link{mapes}."
    # fadmits<-.self$predict_level(n.ahead=n.ahead, sea.on=FALSE)
      sea<-.self$predict_level(n.ahead=n.ahead, sea.on=TRUE)
      
      idx.dates <- (index(Y) >=end.date)
      data_validation<-na.omit(add_daily_ldl(Y[idx.dates], LeadIndCol = LeadIndCol))[1:n.ahead]
      
      newTarg_validation<-data_validation[,c("newTarg")]
      compare<-cbind(newTarg_validation, sea[,1]) #fadmits[,1]
      names(compare)<-c("Actual", "Forecast")
      
      # mape.trend <- 100*(abs(compare$Actual - compare$ForecastTrend)/
      #                      compare$Actual) %>% mean
      mape.sea <- 100*(abs(compare$Actual - compare$Forecast)/compare$Actual) %>%
        mean
      
      mae<-abs(compare$Actual - compare$Forecast) %>% mean
      rmse<-sqrt(mean((compare$Actual - compare$Forecast)^2))
      coverage<-100*sum(and(sea[,2]<=compare$Actual, sea[,3]>=compare$Actual))/n.ahead
      
    return(list(mape=mape.sea, mae=mae, rmse=rmse, coverage=coverage))
  }
)
)
