# Created by: Edwin Tang
# Created on: 07/11/2024

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

#' @title Supply exogenous predictors for new predictions 
#
#' @description Accessor method to provide new exogenous predictors for the 
#' prediction period.
#'
#' @param object FilterResults or FilterReusltsLI object
#' @param new.xts An xts object containing new exogenous predictors
#' @param idx Character string, either "lead" or "target", for which exogenous variables are supplied. 
#' Only applicable for FilterResultsLI object. Defaults to NULL.
#'
#' @importFrom xts is.xts
#' @examples
#' library(tsgc)
#' 
#' #FilterResults example
#' #Load Gauteng weather 
#' data(gauteng_weather_2021, package = "tsgc")
#' gauteng_weather<-gauteng_weather_2021[,c(1,3)]
#' 
#' # Set up model and estimate it
#' model_weather <- SSModelDynamicGompertz$new(Y = gauteng, xpred=gauteng_weather,
#'                                             start.date=as.Date("2021-02-01"), 
#'                                            end.date=as.Date("2021-04-19"))
#' res_weather <- estimate(model_weather)
#' res_weather$xpred.new
#' 
#' # Feed future weather data into the results object. Subsetting of gauteng_weather 
#' #is done inside the function.
#' supply_xpred.new(res_weather,gauteng_weather)
#' res_weather$xpred.new
#' 
#' #FilterResultsLI example
#' xpred_lead<-xpred_targ<-england_weather_2021[,1:4]
#' mod<-SSModelLeadingIndicator$new(england[,1:2], n.lag=4, xpred_lead=xpred_lead, xpred_targ=xpred_targ, 
#'                                 start.date = as.Date("2021-04-30"), 
#'                                 end.date = as.Date("2021-07-24"))
#' res_lead.x<-estimate(mod)
#'
#' supply_xpred.new(res_lead.x,england_weather_2021[,1:4],idx='lead')
#' supply_xpred.new(res_lead.x,england_weather_2021[,1:4],idx='targ')
#' res_lead.x
#' 
#' @export
supply_xpred.new<-function(object, new.xts, idx=NULL){
  if (!is.xts(new.xts)){
    stop("new.xts is not an xts object.")
  }
  if (class(object)=="FilterResultsLI"){
    if (idx=="lead"){
      object$xpred_lead.new<-new.xts
      print("xpred_lead.new registered.")
    } else if (idx=="targ"){
      object$xpred_targ.new<-new.xts
      print("xpred_targ.new registered.")
    } else {
      stop("Please specify idx, which is either 'lead' or 'targ'.")
    }
  } else if (class(object)=="FilterResults"){
    object$xpred.new<-new.xts
    print("xpred.new registered.")
  } else {
    stop("Object is not of class 'FilterResultsLI' or 'FilterResults'.")
  }
}

#' @title Extract output of FilterResults or FilterResultsLI
#' 
#' @description Accessor method to access the fitted KFS model from `FilterResults`
#'
#' @param object FilterResults object
#' 
#' @returns The fitted KFS model
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Return KFS object in output of res
#' output(res)
#'
output<-function(object){
  return(object$output)
}

#' @title Extract SSModel object within KFS object
#
#' @description Accessor method to access the fitted SSModel
#'
#' @param object KFS object
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, end.date=as.Date("2020-07-06"), q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#' 
#' # Extract Z matrix from output(res)
#' modelKFS(output(res))
#'
modelKFS<-function(object){
  return(object$model)
}

#' @title Extract number of seasonal components used in KFS
#
#' @description Accessor method to access the number of seasonal component used in KFS object
#'
#' @param object KFS object
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' #Return number of seasonal components
#' seasonalComp(output(res))
#'
seasonalComp<-function(object){
  attr(modelKFS(object)$terms, "specials")$SSMseasonal
}

#' @title Extract filtered state estimates used in KFS
#
#' @description Accessor method to access filtered state estimate in KFS object
#'
#' @param object KFS object
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' #Return filtered state estimates
#' att(output(res))
#'
#' @returns Filtered state estimate in KFS object
#' 
att<-function(object){
  object$att
}

#' @title Extract error covariance matrix of filtered states from KFS
#
#' @description Accessor method to access the non-diffuse part of the error 
#' covariance matrix of filtered states from KFS object
#'
#' @param object KFS object
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' #Return covariance matrix
#' Ptt(output(res))
#'
Ptt<-function(object){
  object$Ptt
}

#' @title Extract error covariance matrix of smoothed states from KFS
#
#' @description Accessor method to access the non-diffuse part of the error 
#' covariance matrix of smoothed states from KFS object
#'
#' @param object KFS object
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' #Return covariance matrix
#' get_V(output(res))
#'
get_V<-function(object){
  object$V
}

#' @title Extract matrices used in observation, state and disturbance equation 
#' in KFS object
#
#' @description Accessor method to access the matrices used in observation, 
#' state and disturbance equation in KFS object
#'
#' @param object KFS object
#' @param matrix String, indicating a matrix component of SSModel, e.g. H,T,R,Q
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
#' # Extract Z matrix from output(res)
#' matrixKFS(output(res),"Z")
#' 
matrixKFS<-function(object,matrix){
  modelKFS(object)[[matrix]]
}

#' @title Extract time series y in SSModel
#
#' @description Accessor method to access time series y in `SSModel` object
#'
#' @param object An `SSModel` object
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' #Return y
#' gety(modelKFS(output(res)))
#'
gety<-function(object){
  object$y
}

#' @title Extract prediction y.hat in predict.all output
#
#' @description Accessor method to access prediction y.hat in `SSModel` object
#'
#' @param object An object from the output of predict.all
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' 
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, 
#'                                     end.date=as.Date("2020-07-20"))
#' 
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Get object from predict.all
#' all_predictions<-res$predict_all(n.ahead=7)
#' 
#' # Get prediction
#' gety.hat(all_predictions)
#' 
gety.hat<-function(object){
  object$y.hat
}


#' @title Extract alphahat in KFS object
#
#' @description Accessor method to access alphahat (smoothed state estimates) 
#' from a fitted `KFS` object which has had the Kalman filter applied to it
#'
#' @param object A `KFS` object
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' #Return alphahat
#' alphahat(output(res))
#'
alphahat<-function(object){
  object$alphahat
}

#' @title Calling estimate method for SSModelDynamicGompertz or SSModelLeadingIndicator class
#'
#' @description Accessor method to obtain estimated model for `SSModelDynamicGompertz` class
#'
#' @param model A `SSModelDynamicGompertz` or `SSModelLeadingIndicator` object
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' estimate(model)
#' 
#' @export
estimate<-function(model){
   model$estimate()
}

#' @title Calling print method for classes in tsgc
#'
#' @description Accessor method to print a short description for the objects of
#' `SSModelLeadingIndicator` class
#'
#' @param model A `SSModelLeadingIndicator` object
#' 
#' @method print SSModelLeadingIndicator
#' 
#' @examples
#' library(tsgc)
#' 
#' # Specify a model
#' out_eng <- tsgc::SSModelLeadingIndicator(
#' Y = england[, 1:2], n.lag = 4, sea.period = 7,
#' start.date = as.Date("2021-04-30"), end.date = as.Date("2021-07-24"))
#' 
#' # Print a short description of the model object
#' print(out_eng)
#' 
#' 
#' @export
print.SSModelLeadingIndicator <- function(model) {
  # Call the object's print() method if it exists
  if (!is.null(model$print) && is.function(model$print)) {
    model$print()
  } else {
    stop("The object does not have a valid 'print' method.")
  }
}

#' @title Calling plot method for `SSModelLeadingIndicator` class
#'
#' @description Accessor method to call the plot method of an object of 
#' `SSModelLeadingIndicator` class
#'
#' @param model A `SSModelLeadingIndicator` object
#' @method plot SSModelLeadingIndicator
#' 
#' @examples
#' library(tsgc)
#' 
#' # Specify a model
#' out_eng <- tsgc::SSModelLeadingIndicator(
#' Y = england[, 1:2], n.lag = 4, sea.period = 7,
#' start.date = as.Date("2021-04-30"), end.date = as.Date("2021-07-24"))
#' 
#' plot(out_eng)
#' 
#' @export
plot.SSModelLeadingIndicator <- function(model,title=NULL, series.name.lead="Leading Indicator", 
                                         series.name.target="Target Variable",
                                         date_break=NULL, take.log=TRUE) {
  # Call the object's plot() method if it exists
  if (!is.null(model$plot) && is.function(model$plot)) {
    model$plot(title, series.name.lead, 
               series.name.target,date_break,take.log)
  } else {
    stop("The object does not have a valid 'plot' method.")
  }
}


#' @title Calling summary method for classes in tsgc
#'
#' @description Accessor method to show a summary for the objects of
#' `SSModelLeadingIndicator` class
#'
#' @param model A `SSModelLeadingIndicator` object
#' @method summary SSModelLeadingIndicator
#' 
#' @examples
#' library(tsgc)
#' 
#' # Specify a model
#' out_eng <- tsgc::SSModelLeadingIndicator(
#' Y = england[, 1:2], n.lag = 4, sea.period = 7,
#' start.date = as.Date("2021-04-30"), end.date = as.Date("2021-07-24"))
#' 
#' summary(out_eng)
#' 
#' @export
summary.SSModelLeadingIndicator <- function(model) {
  # Call the object's summary() method if it exists
  if (!is.null(model$summary) && is.function(model$summary)) {
    model$summary()
  } else {
    stop("The object does not have a valid 'summary' method.")
  }
}

#' @title Calling print method for SSModelDynamicGompertz class
#'
#' @description Accessor method to print a short description for the objects of
#' `SSModelDynamicGompertz` class
#'
#' @param model A `SSModelDynamicGompertz` object
#' @method print SSModelDynamicGompertz
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' 
#' # Print a short description of the model object
#' print(model)
#' 
#' @export
print.SSModelDynamicGompertz <- function(model) {
  # Call the object's print() method if it exists
  if (!is.null(model$print) && is.function(model$print)) {
    model$print()
  } else {
    stop("The object does not have a valid 'print' method.")
  }
}

#' @title Calling summary method for SSModelDynamicGompertz class
#'
#' @description Accessor method to show a summary for the objects of
#' `SSModelDynamicGompertz` class
#'
#' @param model A `SSModelDynamicGompertz` object
#' @method summary SSModelDynamicGompertz
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
#' @export
summary.SSModelDynamicGompertz <- function(model) {
  # Call the object's summary() method if it exists
  if (!is.null(model$summary) && is.function(model$summary)) {
    model$summary()
  } else {
    stop("The object does not have a valid 'summary' method.")
  }
}



#' @title Calling plot method for SSModelDynamicGompertz class
#'
#' @description Accessor method to call the plot method of an object of 
#' `SSModelDynamicGompertz` class
#'
#' @param model A `SSModelDynamicGompertz` object
#' @method plot SSModelDynamicGompertz
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
#' plot(model)
#' 
#' @export
plot.SSModelDynamicGompertz <- function(model,title=NULL, 
                                        series.name="target variable", 
                                        date_break=NULL, MA_period=7) {
  # Call the object's plot() method if it exists
  if (!is.null(model$plot) && is.function(model$plot)) {
    model$plot(title, series.name, date_break, MA_period)
  } else {
    stop("The object does not have a valid 'plot' method.")
  }
}

#' @title Calling summary method for FilterResults
#'
#' @description Accessor method to show a summary for the objects of
#' `FilterResults` class
#'
#' @param model A `FilterResults` object
#' @method summary FilterResults
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Return KFS object in output of res
#' summary(res)
#' 
#' @export
summary.FilterResults <- function(model) {
  # Call the object's summary() method if it exists
  if (!is.null(model$summary) && is.function(model$summary)) {
    model$summary()
  } else {
    stop("The object does not have a valid 'summary' method.")
  }
}

#' @title Calling print method for FilterResults class
#'
#' @description Accessor method to print a short description for the objects of
#' `FilterResults` class
#'
#' @param model A `FilterResults` object
#' @method print FilterResults
#' 
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng, q = 0.005, end.date=as.Date("2020-07-20"))
#' # Estimate a specified model
#' res <- estimate(model)
#' 
#' # Return short description of fitted model
#' print(res)
#' 
#' @export
print.FilterResults <- function(model) {
  # Call the object's print() method if it exists
  if (!is.null(model$print) && is.function(model$print)) {
    model$print()
  } else {
    stop("The object does not have a valid 'print' method.")
  }
}

#' @title Calling summary method for FilterResultsLI
#'
#' @description Accessor method to show a summary for the objects of
#' `FilterResultsLI` class
#'
#' @param model A `FilterResultsLI` object
#' @method summary FilterResultsLI
#' 
#' @examples
#' library(tsgc)
#' 
#' out_eng <- tsgc::SSModelLeadingIndicator(
#' Y = england[, 1:2], n.lag = 4, sea.period = 7,
#' start.date = as.Date("2021-04-30"), end.date = as.Date("2021-07-24"))
#' 
#' res_eng<-estimate(out_eng)
#' summary(res_eng)
#' 
#' @export
summary.FilterResultsLI <- function(model) {
  # Call the object's summary() method if it exists
  if (!is.null(model$summary) && is.function(model$summary)) {
    model$summary()
  } else {
    stop("The object does not have a valid 'summary' method.")
  }
}

#' @title Calling print method for FilterResultsLI class
#'
#' @description Accessor method to print a short description for the objects of
#' `FilterResultsLI` class
#'
#' @param model A `FilterResultsLI` object
#' @method print FilterResultsLI
#' 
#' @examples
#' library(tsgc)
#' 
#' out_eng <- tsgc::SSModelLeadingIndicator(
#' Y = england[, 1:2], n.lag = 4, sea.period = 7,
#' start.date = as.Date("2021-04-30"), end.date = as.Date("2021-07-24"))
#' 
#' res_eng<-estimate(out_eng)
#' print(res_eng) 
#' 
#' @export
print.FilterResultsLI <- function(model) {
  # Call the object's print() method if it exists
  if (!is.null(model$print) && is.function(model$print)) {
    model$print()
  } else {
    stop("The object does not have a valid 'print' method.")
  }
}

