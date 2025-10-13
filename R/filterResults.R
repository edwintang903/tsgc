setOldClass("KFS")
#'
#' @title Class for estimated Dynamic Gompertz Curve model
#'
#' @description Class for estimated Dynamic Gompertz Curve model and contains
#' methods to extract smoothed/filtered estimates of the states, the level of
#' the incidence variable \eqn{y}, and forecasts of \eqn{y}. The output from the estimate method
#' of the SSModelDynGompertz class is of the class FilterResults.
#' 
#' @field data_xts An xts object containing the non-reinitialized cummulated 
#' variable.
#' @field xpred_logical Logical value indicating whether exogenous predictors were 
#' used to estimate the FilterResults object. 
#' @field index The list of dates in the index of \code{data_xts}.
#' @field reinit.date The reinitialisation date of the estimated \code{SSModelDynamicGompertz} model (if applicable). 
#' Should be specified as an object of class \code{"Date"}.
#' @field ar1 Logical value indicating whether an ar1 component should be 
#' included in the model.
#' @field output A \code{KFS} results object obtained after fitting a 
#' \code{SSModelDynamicGompertz} model.
#' @field xpred.new An xts object containing exogenous predictors to be used in 
#' prediction. Defaults to \code{NULL}, and should be provided if xpred is 
#' used for model estimation.
#' @field resolution A character object showing the time resolution of the data 
#' in \code{data_xts}. Options are "daily", "monthly, "quarterly" and "yearly".
#' Automatically estimated when \code{data_xts} is provided.
#' 
#' @references Harvey, A. C. and Kattuman, P. (2021). A Farewell to R:
#' Time Series Models for Tracking and
#' Forecasting Epidemics, Journal of the Royal Society Interface, vol 18(182):
#' 20210179
#'
#' @importFrom xts periodicity last
#' @importFrom magrittr %>%
#' @importFrom methods new
#' @importFrom abind abind
#' @importFrom zoo as.yearqtr as.yearmon
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#'
#' # Estimation and prediction settings
#' estimation.date.end=as.Date("2020-07-20")
#' plt.length=30
#' 
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, 
#' end.date=estimation.date.end)
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
#' # Forecast 7 days ahead from the end of the estimation window
#' res$predict_level(n.ahead = 7,
#'   confidence.level = 0.68, sea.on=TRUE)
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
#' # Plot forecast and realised log growth rate of cumulative cases
#' res$plot_log_forecast(Y=gauteng,n.ahead=7,
#' plt.start.date=estimation.date.end-plt.length)
#' 
#' # Plot forecast of new cases 7 days ahead
#' res$plot_forecast(n.ahead=7,
#' plt.start.date = estimation.date.end-plt.length,
#' series.name="hospitalizations")
#' 
#' # Plot forecasts and outcomes over evaluation period
#' res$plot_holdout(Y=gauteng,n.ahead=7, series.name="hospitalizations")
#' 
#' # Plot filtered gy, g and gamma
#' res$plot_gy_components()
#' 
#' # Plot filtered gy, g and gamma
#' res$plot_gy_ci()
#' 
#' # Return MAPE of forecast
#' res$mapes(n.ahead=7,gauteng)
#'
#' @export
#'
FilterResults <- setRefClass(
  "FilterResults",
  fields = list(
    data_xts = "xts",
    xpred_logical = "ANY",
    xpred.new="ANY",
    index = "ANY",
    reinit.date= "ANY",
    ar1 = "logical",
    output = "KFS",
    sea.period="numeric",
    resolution="character"),
  methods = list(
    initialize = function(data_xts,xpred_logical,index,reinit.date, ar1, 
                          output, sea.period, xpred.new=NULL, resolution="daily")
    {
      "Create an instance of the \\code{FilterResults} class with fields defined
      earlier in the fields section."
      data_xts<<-data_xts
      index <<- index
      xpred_logical<<-xpred_logical
      xpred.new<<-xpred.new
      reinit.date<<-reinit.date
      ar1<<-ar1
      output <<- output
      sea.period<<-sea.period
      resolution<<-get_time_resolution(index)
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
      \\subsection{Return Value}{\\code{xts} object containing the point
      forecasts and upper and lower bounds of
      the forecast interval.}"
      if (!is.null(reinit.date)){
        y.cum<-reinitialise_dataframe(data_xts, reinit.date)
      } else {
        y.cum<-data_xts
      }
      model <- modelKFS(output)
      n <- attr(model, "n")
      p <- attr(model, "p")

      freq <- unclass(periodicity(y.cum))$label
      endtime <- end(gety(model)) + c(0, n.ahead)
      filtered.out <- .self$predict_all(n.ahead, sea.on = sea.on,
                                        return.all = FALSE, 
                                        confidence.level = confidence.level)

      # # 1. Extract parameters.
      timespan <- n + 0:n.ahead

      # Calculate g.t as exponent of y.t
      g.t <- exp(gety.hat(filtered.out)[,1])
      g.t.lwr <- exp(gety.hat(filtered.out)[,2])
      g.t.upr <- exp(gety.hat(filtered.out)[,3])

      # Forecast dates
      v_dates_end <- if (resolution=='daily'){
        seq(last(index(y.cum)), last(index(gety.hat(filtered.out))), by = freq)
      } else if (resolution=='quarterly'){
        as.yearqtr(seq(as.numeric(last(index(y.cum))),
                                as.numeric(last(index(gety.hat(filtered.out)))),
                                by=0.25))
      } else if (resolution=='yearly'){
        as.yearmon(seq(as.numeric(last(index(y.cum))),
                       as.numeric(last(index(gety.hat(filtered.out)))),
                       by=1))
      } else if (resolution=='monthly'){
        as.yearmon(seq(as.numeric(last(index(y.cum))),
                       as.numeric(last(index(gety.hat(filtered.out)))),
                       by=1/12))}

      y.hat <- xts(matrix(NA, nrow = n.ahead + 1, ncol = 3),
                   order.by = v_dates_end)
      y.hat[v_dates_end[1],] <- y.cum[v_dates_end[1]]
      for (i in seq_len(length(v_dates_end[-1]))) {
        date.forecast <- v_dates_end[i + 1]
        date.lag <- if (resolution=='daily' || resolution=='yearly'){date.forecast - 1} 
        else if (resolution=='quarterly'){date.forecast - 0.25}
        else if (resolution=='monthly'){date.forecast - 1/12}
        # Update level
        y.hat[date.forecast, 1] <- as.numeric(y.hat[date.lag, 1]) *
          as.numeric(1 + g.t[date.forecast,])

        # Make prediction intervals
        y.hat[date.forecast, 2] <- as.numeric(y.hat[date.lag, 1]) *
          as.numeric(1 + g.t.lwr[date.forecast,])
        y.hat[date.forecast, 3] <- as.numeric(y.hat[date.lag, 1]) *
          as.numeric(1 + g.t.upr[date.forecast,])
      }

      # Difference output if requested
      d <- if (return.diff) { diff(y.hat[, 1])[-1] } else { (y.hat[, 1])[-1] }

      ci_bounds <- if (return.diff) {
        (y.hat[, 2:3] - as.vector(y.hat[, 1]))[-1] + as.vector(d)
      } else { y.hat[2:dim(y.hat)[1], 2:3] }

      pred <- vector("list", length = p)
      pred[[p]] <- cbind(fit = d, prediction = ci_bounds)
      pred <- lapply(pred, ts, end = endtime, frequency = 1)

      y.hat <- xts(pred[[p]], order.by = v_dates_end[-1])
      names(y.hat)[2:3] <- list('lower', 'upper')

      return(as.xts(y.hat))
    },
    print_estimation_results = function() {
      "Prints a table of estimated parameters in a format ready to paste into
      LaTeX."
      H <- output$model$H[, , 1]
      Q_gamma <- output$model$Q[2, 2, 1]
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
      \\subsection{Return Value}{\\code{xts} object containing the forecast
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
          xpred.new<<-get_timeframe(xpred.new,tail(index,1)+1,tail(index,1)+n.ahead)
          
          newZ<-array(new.model$Z[,,dim(new.model$Z)[3]], 
                      dim = c(dim(new.model$Z)[1], dim(new.model$Z)[2], n.ahead))
          newZ[,1:dim(xpred.new)[2],]<-t(xpred.new)
          
          new.model$Z <- abind::abind(
            new.model$Z,
            newZ,
            along = 3
          )
          
          model_output <- KFS(new.model)
          new.Q <- new.model$Q
          if (ar1){
            #AR1 with sea.period
            if (!is.null(sea.period)){
              ar1_index<-dim(new.Q)[1]
              newdata<-SSModel(rep(NA,dim(xpred.new)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMseasonal(
                                 period = sea.period, 
                                 Q = new.Q[3,3,1],
                                 sea.type = "trigonometric")
                               +SSMregression(~xpred.new)
                               +SSMcustom(Z=1,T=1,R=1,Q=new.Q[ar1_index,ar1_index,1],state_names="ar1"))
            } else {
              #AR1 and no sea.period
              ar1_index<-dim(new.Q)[1]
              newdata<-SSModel(rep(NA,dim(xpred.new)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMregression(~xpred.new)
                               +SSMcustom(Z=1,T=1,R=1,Q=new.Q[ar1_index,ar1_index,1],state_names="ar1"))
            }
          } else {
            #sea period only
            if (!is.null(sea.period)){
              newdata<-SSModel(rep(NA,dim(xpred.new)[1])
                               ~SSMtrend(degree = 2,
                                         Q = list(matrix(0), matrix(new.Q[2,2,1])))
                               +SSMseasonal(
                                 period = sea.period, 
                                 Q = new.Q[3,3,1],
                                 sea.type = "trigonometric")
                               +SSMregression(~xpred.new))
            } else {
              #no sea period 
            newdata<-SSModel(rep(NA,dim(xpred.new)[1])
                             ~SSMtrend(degree = 2,
                                       Q = list(matrix(0), matrix(new.Q[2,2,1])))
                             +SSMregression(~xpred.new))
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
      
      dates <- seq_dates(from=index[1], length.out = (oldn + n.ahead), 
                         resolution=resolution)

      y.hat <- xts::xts(
        rbind(y.t.t, y.hat.kfas),
        order.by = dates)
      names(y.hat)<-c("y.hat","y.hat.upr","y.hat.lwr")
      
      i.level <- grep("level", colnames(att(model_output)))
      level.t.t <- xts::xts(att(model_output)[, i.level], order.by = dates) %>%
        as.xts()
      names(level.t.t)<-c("level.t.t")
      
      i.slope <- grep("slope", colnames(att(model_output)))
      slope.t.t <- xts::xts(att(model_output)[, i.slope], order.by = dates) %>%
        as.xts()
      names(slope.t.t)<-c("slope.t.t")

      if (!return.all) {
        y.hat <- y.hat %>%
          subset(index(.) > tail(index, 1))
        level.t.t <- level.t.t %>%
          subset(index(.) > tail(index, 1))
        slope.t.t <- slope.t.t %>%
          subset(index(.) > tail(index, 1))
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
      \\subsection{Return Value}{\\code{xts} object containing
      smoothed/filtered growth rates and components (\\eqn{\\delta} and
      \\eqn{\\gamma}), where applicable.}"
      kfs_out <- output
      idx <- index

      if (smoothed) {
        att <- alphahat(kfs_out)
      } else {
        att <- att(kfs_out)
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

      kfs_out <- output
      idx <- index

      if (smoothed) {
        att <- alphahat(kfs_out)
        var <- get_V(kfs_out)
      } else {
        att <- att(kfs_out)
        var <- Ptt(kfs_out)
      }

      filtered_slope <- xts(att[, "slope"], order.by = idx)
      filtered.level <- xts(att[, "level"], order.by = idx)
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope

      idx.slope <- grep("slope", colnames(att(kfs_out)))
      ci <- qnorm((1 - confidence.level) / 2) *
        sqrt(var[idx.slope, idx.slope,]) %o% c(1, -1)
      ci_bounds <- as.vector(gy.t) + ci

      pred <- xts(cbind(gy.t, ci_bounds), order.by = idx)
      colnames(pred) <- c("fit","lower","upper")

      return(pred)
    },
    print=function(){
      "Provides a quick glimpse of model states and standard errors."
      cat("Object of FilterResults Class\n")
      cat("  - Model States and Standard Errors\n")
      base::print(output)
    },
    summary=function(){
      "Supplies details of the filterResults object, such as estimated 
      parameter values, start and end dates of estimation."
      H <- matrixKFS(output, "H")[, , 1]
      Q_gamma <- matrixKFS(output, "Q")[2, 2, 1]
      if (sea.period>1){  
        Q_seasonal <- matrixKFS(output, "Q")[3, 3, 1]
      }
      
      start.date <- index[1]
      end.date <- index[length(index)]
      
      cat("Summary of FilterResults Object\n")
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
    plot_forecast=function(n.ahead=14, confidence.level = 0.68, 
                            title=NULL, plt.start.date=NULL, 
                            series.name="target variable") {
      "Generates a forecast plot for the difference in the cumulative variable,
      showing actual values, forecasts including seasonal components,
      and prediction intervals around the forecasts. 
      For more details, see \\link{plot_forecast}."
      
      if (xpred_logical){
        if (is.null(xpred.new)){
          stop("xpred.new cannot be NULL.")
        } 
      }
      
      Date <- Data <- Forecast <- ForecastTrend <- lower <- upper <- NULL
      if (is.null(title)) {title <- ""}
       
      estimation.date.end <- tail(index, 1)
      
      if (!is.null(reinit.date)){
        Y<-reinitialise_dataframe(data_xts, reinit.date)
      } else{
        Y<-data_xts
      }
    y.level.est <- Y[index]
    if (is.null(plt.start.date)) {plt.start.date <- head(index, 1)}
    
    y.hat.diff.final.ci <- .self$predict_level(
      n.ahead = n.ahead, confidence.level = confidence.level,
      sea.on=TRUE
    )
    # y.hat.diff.final <- .self$predict_level(
    #   n.ahead = n.ahead, confidence.level = confidence.level,
    #   sea.on = TRUE
    # )
    # 
    tmp.date <- if (resolution=='daily'){
      as.Date(plt.start.date)
    } else if (resolution=='quarterly' || resolution=='monthly' || resolution=='yearly'){
      as.Date(format(as.yearmon(plt.start.date), format="%Y-%m-%d"))
    }
    s <- sprintf("%s/", format(tmp.date, "%Y-%m-%d"))
    d.plot <- cbind(
      diff(y.level.est)[s],
      #y.hat.diff.final[, 1],
      y.hat.diff.final.ci[, 1]
    )
    names(d.plot) <- c('Data', 'Forecast')
    
    date_col<-if(resolution=='daily'){
      as.Date(index(y.hat.diff.final.ci))} 
    else if (resolution=='quarterly' || resolution=='monthly' || resolution=='yearly') {
      qtr2date(index(y.hat.diff.final.ci))
      }
    
    ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                              date_col))
    colnames(ci) <- c('lower', 'upper', 'date')
    ci[, 'date'] <- as.Date(ci[, 'date'], origin = "1970-01-01")
    
    df_plot <- as.data.frame(d.plot)
    
    df_plot$Date<-if (resolution=='quarterly'){
      qtr2date(as.yearqtr(rownames(df_plot)))
    } else if (resolution=='monthly'|| resolution=='yearly'){
      qtr2date(as.yearmon(rownames(df_plot)))
    } else {
      as.Date(rownames(df_plot))
    }
    
    ggplot2::ggplot(data = df_plot, aes(x = Date)) +
      ggplot2::geom_line(aes(y = Data, color = "Data"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
      ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
      ggplot2::geom_ribbon(data = ci, aes(x = date, ymin = lower, ymax = upper),
                           linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
      labs(x = "Date", y = paste("New", series.name), title = title) +
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
    },
    plot_log_forecast=function(Y,n.ahead = 14, plt.start.date=NULL, title="", caption = "") {
      "Plots actual and filtered values of the log cumulative growth rate 
      (\\eqn{\\ln(g_t)}) in the estimation sample and the forecast and realised 
      log cumulative growth rate out of the estimation sample. For more details,
      see \\link{plot_log_forecast}."
      model <- modelKFS(output)
      if (!is.null(reinit.date)){
        y.eval <- Y %>%
          reinitialise_dataframe(., reinit.date) %>%
          df2ldl() %>%
          subset(index(.) > tail(index,1))
      } else {y.eval<-Y %>%
        df2ldl()%>% subset(zoo::index(.) > tail(index,1))}
      
      y <- xts::xts(model$y %>% as.numeric(), order.by = index)
      p <- attr(model, 'p')
      
      firstpred<-if (resolution=='quarterly'){
        tail(index,1)+0.25
      } else if (resolution=='monthly'){
        tail(index,1)+1/12
      } else {
        tail(index,1)+1
      }
      
      y.hat.all <- .self$predict_all(n.ahead, return.all = TRUE)
      y.pred <-  get_timeframe(y.hat.all$y.hat$y.hat, firstpred)
      filtered.level <- y.hat.all$level
      
      if (p == 1) {
        EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
        
        if (xpred_logical){
          d <- cbind(y, y.pred, get_timeframe(y.eval, firstpred))
          if (!is.null(plt.start.date)) { d <- d[index(d) > plt.start.date] }
          d <- d[index(d) <= tail(index(y.pred),1)]
          names(d) <- c('EstimationSample', 'Forecast', 'RealisedData')
        } else {
          d <- cbind(y, filtered.level, y.pred, get_timeframe(y.eval, firstpred))
          if (!is.null(plt.start.date)) { d <- d[index(d) > plt.start.date] }
          d <- d[index(d) <= tail(index(y.pred),1)]
          names(d) <- c('EstimationSample', 'FilteredLevel', 'Forecast', 'RealisedData')
        }
        
        df_plot <- as.data.frame(d)
        
        df_plot$Date<-if (resolution=='quarterly'){
          qtr2date(as.yearqtr(rownames(df_plot)))
        } else if (resolution=='monthly' || resolution=="yearly"){
          qtr2date(as.yearmon(rownames(df_plot)))
        } else {
          as.Date(rownames(df_plot))
        }
        
        
        if (!xpred_logical){
          color_values <- c("Estimation\nSample" = 1, "Filtered\nLevel" = 2, 
                            "Forecast" = 3, "Realised\nData" = "grey")
          linetype_values <-c("solid",
                              "solid",
                              "solid",
                              "dashed")
          p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date))+
            ggplot2::geom_line(aes(
              y = EstimationSample, color = "Estimation\nSample"), lwd = 0.85) +
            ggplot2::geom_line(aes(y = FilteredLevel, color = "Filtered\nLevel"), lwd = 0.85)
        } else {
          color_values <-c("Estimation\nSample" = 1,
            "Forecast" = 3, "Realised\nData" = "grey")
          linetype_values<-c("solid", "solid", "dashed")
          p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date))+
            ggplot2::geom_line(aes(
              y = EstimationSample, color = "Estimation\nSample"), lwd = 0.85)
        }
        p1 <- p1+
          ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
          ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData"),
                             lwd = 0.85) +
          ggplot2::scale_color_manual(values = color_values) +
          scale_linetype_manual(values = linetype_values) +
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
            plot.caption = element_text(size = rel(1))
          )
      } else if (p == 2) {
        g_1 <- g_2 <- delta <- Forecast <- RealisedData <- NULL
        d <- cbind(y, filtered.level, y.pred[,2],
                   y.eval[index(y.eval)>tail(index(index),1),2])
        d <- d[index(d) <= tail(index(y.pred),1)]
        names(d) <- c('g_1', 'g_2', 'delta', 'Forecast', 'RealisedData')
        
        df_plot <- as.data.frame(d)
        df_plot$Date <- if (resolution=='quarterly'){
          qtr2date(as.yearqtr(rownames(df_plot)))
        } else if (resolution=='monthly' || resolution=="yearly"){
          qtr2date(as.yearmon(rownames(df_plot)))
        } else {
          as.Date(rownames(df_plot))
        }
        
        p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
          ggplot2::geom_line(aes(y = g_1, color = "g_1")) +
          ggplot2::geom_line(aes(y = g_2, color = "g_2")) +
          ggplot2::geom_line(aes(y = g_2, color = "delta")) +
          ggplot2::geom_line(aes(y = Forecast, color = "Forecast")) +
          ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData")) +
          ggplot2::scale_color_manual(
            values = c(1, 2, 3, 4, 'grey')) +
          ggplot2::scale_linetype_manual(
            values = c("solid", "solid", "solid", "solid", "dashed")
          ) +
          ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
          labs(x = "Date", y = "Log Growth Rate", caption = caption,
               title = title
          ) +
          theme_economist_white(gray_bg = FALSE) +
          scale_fill_economist() +
          theme(legend.title = element_blank()) +
          theme(
            text = element_text(size = rel(1.)),
            axis.text = element_text(size = rel(1)),
            axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
            axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
            plot.title = element_text(margin=margin(b=5)),
            plot.caption = element_text(size = rel(1))
          )
      } else { stop('NotImplemented Error') }
      
      return(p1)
    }, 
    plot_gy_components = function(plt.start.date = NULL,
                                   smoothed = FALSE, title = NULL){
      "Plots the growth rates and slope of the log cumulative growth rate 
      against the dates in estimation sample. 
      For more details, please see \\link{plot_gy_components}."
      Value <- Variable <- NULL
      # Determine plot start date
      if(is.null(plt.start.date)) plt.start.date <-.self$index[1]
      
      # Get gy.t, g.t and gamma
      gy.components <-.self$get_growth_y(return.components = TRUE, smoothed =
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
    plot_gy_ci = function(plt.start.date = NULL, smoothed = FALSE,
                           title = NULL, series.name = NULL, pad.right = NULL){
      "Plots the growth rates and the slope of the log cumulative growth rate 
      against the dates in estimation sample. 
      For more details, please see \\link{plot_gy_ci}."
      Date <- fit <- upper <- lower <- NULL
      
      # Determine plot start date
      if(is.null(plt.start.date)) plt.start.date <-.self$index[1]
      
      # Get confidence intervals to plot
      gy.ci<-.self$get_gy_ci(smoothed = smoothed)
      
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
        end.date <- tail(index(gy.ci),1)
        p1 <- p1 +
          ggplot2::scale_x_date(
            limits = c(as.Date(plt.start.date), end.date + pad.right))
      }
      
      return(p1)
    }, 
    plot_holdout = function(Y, n.ahead=14,
                            confidence.level = 0.68,
                            series.name = "target variable",
                             title= NULL, caption = NULL) {
      "Plots the forecast of new cases (the difference of the cumulated
      variable) over a holdout sample. For more details, please refer to 
      \\link{plot_holdout}."
      
      if (xpred_logical){
        if (is.null(xpred.new)){
          stop("xpred.new cannot be NULL.")
        } 
      }
      
      if (!is.null(reinit.date)){
        Y.est<-reinitialise_dataframe(data_xts, reinit.date)
      } else {
        Y.est<-data_xts
      }
      
      model <- modelKFS(output)
      
      y.level.est <- Y.est[index]
      estimation.date.end <- tail(index, 1)
      
      p <- attr(model, 'p')
      if(p!=1) { stop('NotImplementedError') }
      
      #Evaluation values
      y.eval.diff <- diff(Y) %>% na.omit
      
      y.hat.diff.final.ci <-.self$predict_level(
        n.ahead = n.ahead,  sea.on=TRUE,
        confidence.level = confidence.level
      )
      
      if (resolution=='daily' || resolution=='yearly'){
        ids=(index(y.eval.diff)>estimation.date.end) & 
          (index(y.eval.diff)<estimation.date.end+n.ahead+1)
      } else if (resolution=='quarterly'){
        ids=(index(y.eval.diff)>estimation.date.end) & 
          (index(y.eval.diff)<estimation.date.end+(n.ahead+1)/4)
      } else if (resolution=='monthly'){
        ids=(index(y.eval.diff)>estimation.date.end) & 
          (index(y.eval.diff)<estimation.date.end+(n.ahead+1)/12)
      }
      
      d <- cbind(
        y.eval.diff[ids,],
        y.hat.diff.final.ci[, 1]
      )
      names(d) <- c('Actual', 'Forecast')
      d.eval <- na.omit(d)
      
      df_plot <- as.data.frame(d)
      df_plot$Date <- if (resolution=='quarterly'){
        qtr2date(as.yearqtr(rownames(df_plot)))
      } else if (resolution=='monthly' || resolution=="yearly"){
        qtr2date(as.yearmon(rownames(df_plot)))
      } else {
        as.Date(rownames(df_plot))
      }
      
      if (any(d.eval$Actual==0)){
        warning("Validation data contains zeros. MAPE is not a reliable measure.")
      }
      
      mape.sea <- 100*(abs(d.eval$Actual - d.eval$Forecast)/d.eval$Actual) %>%
        mean %>% round(2)
      smape<-mean(100*(abs(d.eval$Actual - d.eval$Forecast)/(d.eval$Actual+d.eval$Forecast))) %>% round(2)
      mae<-abs(d.eval$Actual - d.eval$Forecast) %>% mean %>% signif(digits=3)
      rmse<-sqrt(mean((d.eval$Actual - d.eval$Forecast)^2)) %>% signif(digits=3)
      
      date_col<-if(resolution=='daily'){
        as.Date(index(y.hat.diff.final.ci))} 
      else if (resolution=='quarterly' || resolution=='monthly' || resolution=='yearly') {
        qtr2date(index(y.hat.diff.final.ci))
      }
      
      ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                                date_col))
      colnames(ci) <- c('lower', 'upper', 'date')
      ci[, 'date'] <- as.Date(ci[, 'date'], origin = "1970-01-01")
      
      p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
        ggplot2::geom_line(aes(y = Actual, color = "Actual"),lwd = 0.85) +
        ggplot2::geom_line(aes(y = Forecast, color = "Forecast"),lwd = 0.85) +
        ggplot2::scale_color_manual(values = c("black", "#AA2045")) +
        ggplot2::geom_ribbon(data = ci, aes(x = date, ymin = lower, ymax = upper),
                             linetype = 0, linewidth = 0, fill = "#AA2045",
                             alpha = 0.1) +
        labs(x = "Date", y = paste("New",series.name), title = title,
             subtitle = paste("MAPE: ",mape.sea,"%; SMAPE: ",smape,"%; MAE: ", mae,"; RMSE: ", rmse,".", sep="")) +
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
        
        estimation.date.end <- tail(index, 1)
        
        y.eval.diff <-diff(Y[seq_dates(estimation.date.end, resolution, length.out=n.ahead+1)]) %>% na.omit
        
        y.hat.diff.final <- .self$predict_level(
          n.ahead = n.ahead, confidence.level =0.68,
          sea.on = TRUE
        )
        
        # Extract the relevant columns
        filtered_y_eval_diff <- y.eval.diff[index(y.eval.diff) > estimation.date.end]
        forecast_column <- y.hat.diff.final[, 1]
        
        #Form dataframe
        df_plot <- data.frame(
          Actual = coredata(filtered_y_eval_diff),  # Extract data from zoo
          Forecast = forecast_column,
          row.names = index(filtered_y_eval_diff)  # Use index as row names
        )
        
        d.eval <- na.omit(df_plot)
        colnames(d.eval)<-c('Actual', 'Forecast')
        
        if (any(d.eval$Actual==0)){
          warning("Validation data contains zeros. MAPE is not a reliable measure.")
        }
        
        mape.sea <- mean(100*(abs(d.eval$Actual - d.eval$Forecast)/d.eval$Actual))
        smape<-mean(100*(abs(d.eval$Actual - d.eval$Forecast)/(d.eval$Actual+d.eval$Forecast)))
        mae<-abs(d.eval$Actual - d.eval$Forecast) %>% mean
        rmse<-sqrt(mean((d.eval$Actual - d.eval$Forecast)^2))
        coverage<-100*sum(and(y.hat.diff.final[,2]<=y.eval.diff, y.hat.diff.final[,3]>=y.eval.diff))/n.ahead
        
        return(list(mape=mape.sea, smape=smape, mae=mae, rmse=rmse, coverage=coverage))
      }
  )
)
