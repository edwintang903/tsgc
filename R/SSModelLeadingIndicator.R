setOldClass("KFS")
setOldClass("idx_series")

#' @title Class for designing a Leading Indicator Model
#'
#' @description A class for specifying the parameters of a leading indicator model. The model 
#' settings are stored in the fields of this object, and the class contains 
#' methods to obtain a FilterResultsLI object for further analysis.
#'
#' @field Y A cumulated \code{idx_series} with 2 columns: a leading indicator and a target
#' variable. Both the target variable and the lagged leading indicator must be strictly
#' increasing within the estimation window.
#' @field q The signal-to-noise ratio (ratio of slope error variance to target variable observation error variance). 
#' Defaults to \code{'NULL'}, in which case no
#'   signal-to-noise ratio will be imposed. Instead, it will be estimated.
#' @field sea.period A positive integer specifying the period of seasonality used in the
#'   trigonometric seasonal component of the model. For example, use \code{7} for daily 
#'   data to model day-of-the-week effects. A value of \code{0} disables the seasonal 
#'   component entirely. The default is \code{7}, which is suitable for capturing 
#'   weekly seasonality in daily time series.
#' @field n.lag Number of integer positions to lag the leading indicator by.
#' @field xpred_lead An \code{idx_series} object containing the values of exogenous
#' variables for the leading indicator. Dataset must contain values for all positions in the 
#' estimation time frame. Defaults to NULL, indicating no exogenous variables are needed 
#' for the leading indicator.
#' @field xpred_targ An \code{idx_series} object containing the values of exogenous
#' variables for the target variable. Dataset must contain values for all positions in the 
#' estimation time frame. Defaults to NULL, indicating no exogenous variables are 
#' needed for the target variable.
#' @field LeadIndCol The column in \code{Y} that contains the leading indicator. 
#' Defaults to 1.
#' @field start Integer position marking the start of the estimation period.
#' @field end Integer position marking the end of the estimation period.
#' @field calendar An optional \code{idx_calendar} object describing how the
#' integer positions used by \code{Y} map back to calendar time. Since the
#' leading indicator and target series in \code{Y} are already combined
#' onto one shared position axis (see \code{Y}), a single calendar
#' correctly applies to both. Purely cosmetic: never consulted by
#' estimation/filtering, only carried along so that plotting can later
#' translate positions to dates. Defaults to \code{NULL}, in which case
#' plots will be labelled with raw integer positions instead of dates. It
#' is the user's responsibility to ensure the leading indicator and target
#' series were aligned onto a shared, consistent calendar before being
#' combined into \code{Y}; this class does not attempt to reconcile two
#' different calendars.
#'
#' @importFrom methods new setRefClass setOldClass
#' @importFrom KFAS SSMtrend SSMseasonal SSModel SSMregression KFS fitSSM
#' @importFrom purrr partial
#'
#' @examples
#' library(tsgc)
#' set.seed(1)
#' lead <- cumsum(rpois(120, 6)) + 1
#' targ <- cumsum(rpois(120, 8)) + 1
#' Y <- idx_series(cbind(lead, targ), start = 1)
#'
#' # Specify a model with the estimation timeframe
#' model <- SSModelLeadingIndicator(Y = Y, n.lag = 5, sea.period = 0,
#'   LeadIndCol = 1, start = 1, end = 100)
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
#' @export SSModelLeadingIndicator
#' @exportClass SSModelLeadingIndicator
SSModelLeadingIndicator <- setRefClass(
  "SSModelLeadingIndicator",
  fields = list(
    Y = "idx_series",
    q = "ANY",
    sea.period= "numeric",
    n.lag = "numeric",
    LeadIndCol ="numeric",
    xpred_lead = "ANY",  
    xpred_targ = "ANY",
    start = "ANY",
    end = "ANY",
    calendar = "ANY"),
  methods = list(
    initialize = function(Y, n.lag, sea.period=7, q = NULL,
                          LeadIndCol=1, xpred_lead=NULL, xpred_targ=NULL,
                          start=idx_range(Y)[1], end=idx_range(Y)[2],
                          calendar=NULL)
    {"Create an instance of the \\code{SSModelLeadingIndicator} class with the 
      fields laid out at the beginning of the documentation."
      if (length(sea.period) != 1 || 
          !isTRUE(all.equal(sea.period, as.integer(sea.period)))||
          sea.period==1 || sea.period<0){
        stop("sea.period must be a non-negative integer that is not 1.")
      } 
      if (!is.null(xpred_lead) && !is_idx_series(xpred_lead)){
        stop("xpred_lead must be NULL or an idx_series object.")
      } 
      if (!is.null(xpred_targ) && !is_idx_series(xpred_targ)){
        stop("xpred_targ must be NULL or an idx_series object.")
      } 
      if (length(LeadIndCol) != 1 || !(LeadIndCol %in% c(1, 2))){
        stop("LeadIndCol must take values 1 or 2.")
      }
      if (!is.null(calendar) && !is_idx_calendar(calendar)){
        stop("calendar must be NULL or an idx_calendar object.")
      }
      Y <<- Y
      q <<- q
      sea.period<<-sea.period
      n.lag <<- n.lag
      LeadIndCol <<- LeadIndCol
      xpred_lead<<-xpred_lead
      xpred_targ<<-xpred_targ
      start<<-start
      end<<-end
      calendar<<-calendar
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
      
      y$newLead <- idx_lag(y$newLead, n.lag)
      y$LDLlead <- idx_lag(y$LDLlead, n.lag)
      y$cLead <- idx_lag(y$cLead, n.lag)
      
      # Combine into a single, position-aligned idx_series (inner join on
      # position, i.e. keep only positions present in every component
      # series - cLead/cTarg span the full range but newLead/LDLlead etc.
      # are one position shorter due to differencing, and now additionally
      # shifted forward by n.lag positions).
      common_pos <- Reduce(intersect, lapply(y, idx_positions))
      if (length(common_pos) == 0) {
        stop("No overlapping positions remain across the leading indicator ",
             "and target series after lagging by n.lag = ", n.lag,
             "; n.lag may be too large relative to the available data.")
      }
      combined_mat <- do.call(cbind, lapply(y, function(s) idx_values(s[common_pos])))
      colnames(combined_mat) <- names(y)
      y_combined <- idx_series(combined_mat, start = common_pos[1])
      
      # Treat +/-Inf (e.g. from log(0) in df2ldl when an increment happens
      # to be recorded as exactly matching the prior level) as missing,
      # then drop any position with a missing value in any column - the
      # idx_series analogue of xts's na.omit().
      finite_rows <- apply(idx_values(y_combined), 1, function(row) all(is.finite(row)))
      keep_pos <- idx_positions(y_combined)[finite_rows]
      if (length(keep_pos) == 0) {
        stop("No positions remain after removing missing/infinite values (check n.lag and the estimation range).")
      }
      # keep_pos may not be contiguous if interior rows were dropped; the
      # idx_series class only supports contiguous ranges, so we require
      # contiguity here (matching the implicit assumption in the original
      # xts-based code, which relied on na.omit() typically only trimming
      # from the ends when only the lag-induced leading NAs are present).
      if (!identical(keep_pos, seq.int(keep_pos[1], keep_pos[length(keep_pos)]))) {
        stop("Missing/infinite values leave gaps in the middle of the series after filtering, which idx_series cannot represent. Consider a different n.lag or estimation range.")
      }
      y_clean <- y_combined[keep_pos]
      
      y.full <- get_timeframe(y_clean, start)
      y.estimate <- get_timeframe(y_clean, start, end)
      
      newLead_col <- idx_values(y.full)[, "newLead"]
      newTarg_col <- idx_values(y.full)[, "newTarg"]
      if (any(newLead_col<=0) || any(newTarg_col<=0)){
        stop("Y must be strictly increasing within the selected timeframe 
        after lagging the leading indicator. If the cumulative 
           values exhibit plateaus it is necessary to add small increments to 
           eliminate flat segments and allow model estimation. This can be done 
           by ensuring the non-cumulated series is strictly positive.")}
      
      est_pos <- idx_positions(y.estimate)
      
      if (!is.null(xpred_lead)){
        xpred_lead <<- idx_lag(xpred_lead, n.lag)
        est_pos <- intersect(est_pos, idx_positions(xpred_lead))
        if (length(est_pos) == 0){
          stop("xpred_lead (after lagging by n.lag) does not overlap with the estimation timeframe.")
        }
      }
      if (!is.null(xpred_targ)){
        est_pos <- intersect(est_pos, idx_positions(xpred_targ))
        if (length(est_pos) == 0){
          stop("xpred_targ does not overlap with the estimation timeframe.")
        }
      }
      if (!identical(est_pos, seq.int(est_pos[1], est_pos[length(est_pos)]))){
        stop("Restricting the estimation timeframe to the available xpred_lead/xpred_targ data leaves gaps, which idx_series cannot represent.")
      }
      y.estimate <- get_timeframe(y.estimate, est_pos[1], tail(est_pos,1))
      data_mat <- idx_values(y.estimate)[, c("LDLlead","LDLtarg"), drop = FALSE]
      # A single (or zero) usable row after lagging/trimming leaves too
      # little information to estimate this two-equation state-space
      # model, and previously surfaced many calls downstream as an
      # opaque SSModel() "Misspecified H" error rather than pointing at
      # the real cause (an estimation window too narrow relative to
      # n.lag, or too little overlap with xpred_lead/xpred_targ).
      if (nrow(data_mat) < 2) {
        stop("SSModelLeadingIndicator: only ", nrow(data_mat), " usable row(s) ",
             "remain in the estimation window after lagging by n.lag = ", n.lag,
             " and trimming missing/non-overlapping positions; at least 2 are ",
             "required to estimate this model. Widen the estimation window ",
             "(start/end) or reduce n.lag.")
      }
      
      if (!is.null(xpred_lead)){
        xpred_lead <<- get_timeframe(xpred_lead, est_pos[1], tail(est_pos,1))
      }
      if (!is.null(xpred_targ)){
        xpred_targ <<- get_timeframe(xpred_targ, est_pos[1], tail(est_pos,1))
      }
      xreg_lead <- if (!is.null(xpred_lead)) idx_values(xpred_lead) else NULL
      xreg_targ <- if (!is.null(xpred_targ)) idx_values(xpred_targ) else NULL
      
      # Update function allowing the signal-to-noise ratio to be targeted.
      # The signal-to-noise ratio is the variance of the trend component of
      # order 'order' (1 = level, 2 = slope, etc) relative to the variance
      # of the irregular of series 'index' (1 = first column, 2 = second).
      updatesn=function(pars, model, snr, order, index){
        if(any(is.na(model$Q))){
          Q <- as.matrix(model$Q[,,1])
          naQd  <- which(is.na(diag(Q)))
          naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))
          Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
          
          diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
          Q[naQd,naQd][naQnd] <- pars[(length(naQd)+1):(length(naQd)+length(naQnd))]
          model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
        }
        if(!identical(model$H,'Omitted') && any(is.na(model$H))){
          H<-as.matrix(model$H[,,1])
          naHd  <- which(is.na(diag(H)))
          naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
          H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0
          diag(H)[naHd] <-
            exp(0.5 * pars[length(naQd)+length(naQnd)+seq_len(length(naHd))])
          H[naHd,naHd][naHnd] <-
            pars[length(naQd)+length(naQnd)+length(naHd)+seq_len(length(naHnd))]
          model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
          model$Q[order,order,1] <- snr*crossprod(H[index,index])
        }
        model
      }
      # Create the state-space model: a common trend and slope (common
      # trend of degree 2), an extra trend [random walk] in LDLtarg only
      # [degree = 1], and a trigonometric seasonal (period = sea.period,
      # if > 1).
      if (sea.period<2){
        if (is.null(xreg_lead)){
          if (is.null(xreg_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xreg_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        } else {
          if (is.null(xreg_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xreg_lead, type="distinct", index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xreg_lead, type="distinct", index=1)+
                             SSMregression(~xreg_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        }
      }
      else {
        if (is.null(xreg_lead)){
          if (is.null(xreg_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xreg_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        } else {
          if (is.null(xreg_targ)){
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xreg_lead, type="distinct", index=1),
                           H = matrix(c(NA,0,0,NA),2,2))
          } else {
            mod <- SSModel(data_mat ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,NA),2,2),type = 'common')+
                             SSMseasonal(sea.period,Q = matrix(c(0,0,0,0),2,2), sea.type='trigonometric', type='distinct')+
                             SSMtrend(degree = 1, Q = matrix(NA),index=1)+
                             SSMregression(~xreg_lead, type="distinct", index=1)+
                             SSMregression(~xreg_targ, type="distinct", index=2),
                           H = matrix(c(NA,0,0,NA),2,2))
          }
        }
      }
      
      # Number of parameters to estimate: the number of NAs in Q and H.
      npar = sum(is.na(mod$Q)) + sum(is.na(mod$H))
      
      if (is.null(q)){
        fit = fitSSM(mod, rep(0,npar))
      }
      else{
        update = updatesn %>% purrr::partial(snr=q,order=2,index=2)
        fit = fitSSM(mod, rep(0,npar), updatefn = update)
      }
      
      out = KFS(fit$model)
      
      results <- FilterResultsLI$new(
        data = y.full,
        output = out,
        n.lag=n.lag,
        sea.period=sea.period,
        LeadIndCol=LeadIndCol,
        xpred_logical=c(!is.null(xpred_lead),!is.null(xpred_targ)),
        start=est_pos[1],
        end=tail(est_pos,1),
        calendar=calendar)
      return(results)},
    summary = function() {
      "Supplies details of the SSModelLeadingIndicator object, such as estimated 
      parameter values, start and end positions of estimation."
      result<-.self$estimate()
      out <- output(result)
      start_pos<-result$start
      end_pos<-result$end
      
      cat("Summary of SSModelLeadingIndicator Model")
      cat("\n")
      cat("--------------------------------------\n")
      cat("Cumulated Variable:\n")
      base::print(head(idx_values(.self$Y)))
      cat("Model Details:\n")
      cat("  - Model Type: Leading Indicator Model")
      cat("\n")
      cat("  - Seasonal Component: ", ifelse(sea.period>1, "Trigonometric", "None"), "\n")
      cat("  - Period of Seasonality: ", ifelse(sea.period>1, sea.period, "N/A"), "\n")
      cat("  - Estimation start position:", start_pos)
      cat("\n")
      cat("  - Estimation end position:", end_pos)
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(out)
    },
    print = function() {
      "Provides a quick description of the SSModelLeadingIndicator object, providing 
      model states and standard errors."
      
      out <- output(.self$estimate()) #KFS object
      cat("SSModelLeadingIndicator Model")
      cat("\n")
      cat("\n")
      cat("Cumulated Variable:\n")
      base::print(head(idx_values(.self$Y)))
      cat("Number of observations:", length(.self$Y))
      cat("\n")
      cat("Seasonal components?",
          ifelse(is.null(seasonalComp(out)),
                 "No","Yes"),"\n")
    }
  )
)