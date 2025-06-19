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

#' @title Extract output of FilterResults or FilterResultsLI
#
#' @description Accessor method to provide new exogenous predictors for the 
#' estimation period
#'
#' @param object FilterResults or FilterReusltsLI object
#' @param new.xts An xts object containing new exogenous predictors
#' @param idx The series number (must be integers 1 or 2) for which exogenous variables are supplied. 
#' Only applicable for FilterResultsLI object. Defaults to NULL.
#'
#' @export
supply_xpred.new<-function(object, new.xts, idx=NULL){
  if (class(object)=="FilterResultsLI"){
    if (idx==1){
      object$xpred1.new<-new.xts
      print("xpred1.new registered.")
    } else if (idx==2){
      object$xpred2.new<-new.xts
      print("xpred2.new registered.")
    } else {
      stop("Please specify idx, which is either 1 or 2.")
    }
  } else if (class(object)=="FilterResults"){
    object$xpred.new<-new.xts
    print("xpred.new registered.")
  } else {
    stop("Object is not of class 'FilterResultsLI' or 'FilterResults'.")
  }
}

#' @description Accessor method to access the fitted KFS model from `FilterResults`
#'
#' @param object FilterResults object
#' 
#' @returns The fitted KFS model
#'
#' @export
output<-function(object){
  return(object$output)
}

#' @title Extract SSModel object within KFS object
#
#' @description Accessor method to access the fitted SSModel
#'
#' @param object KFS object
#'
#'
#' @export
modelKFS<-function(object){
  return(object$model)
}

#' @title Extract number of seasonal components used in KFS
#
#' @description Accessor method to access the number of seasonal component used in KFS object
#'
#' @param object KFS object
#'
#'
#' @export
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
#' @export
att<-function(object){
  object$att
}

#' @title Extract error covariance matrix used in KFS
#
#' @description Accessor method to access the non-diffuse part of the error 
#' covariance matrix in KFS object
#'
#' @param object KFS object
#'
#'
#' @export
Ptt<-function(object){
  object$Ptt
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
#' @export
matrixKFS<-function(object,matrix){
  modelKFS(object)[[matrix]]
}

#' @title Extract time series y in SSModel
#
#' @description Accessor method to access time series y in `SSModel` object
#'
#' @param object An `SSModel` object
#'
#' @export
gety<-function(object){
  object$y
}

#' @title Extract prediction y.hat in SSModel
#
#' @description Accessor method to access prediction y.hat in `SSModel` object
#'
#' @param object An `SSModel` object
#'
#' @export
gety.hat<-function(object){
  object$y.hat
}


#' @title Extract alphahat in SSModel
#
#' @description Accessor method to access alphahat (smoothed state estimates) 
#' from a fitted `SSModel` object which has had the Kalman filter applied to it
#'
#' @param object A `SSModel` object
#'
#' @export
alphahat<-function(object){
  object$alphahat
}

#' @title Calling estimate method for SSModelDynamicGompertz class
#'
#' @description Accessor method to obtain estimated model for `SSModelDynamicGompertz` class
#'
#' @param model A `SSModelDynamicGompertz` or `SSModelLeadingIndicator` object
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
#' @export
print.SSModelLeadingIndicator <- function(model) {
  # Call the object's print() method if it exists
  if (!is.null(model$print) && is.function(model$print)) {
    model$print()
  } else {
    stop("The object does not have a valid 'print' method.")
  }
}

#' @title Calling plot method for classes in tsgc
#'
#' @description Accessor method to call the plot method of an object of 
#' `SSModelLeadingIndicator` class
#'
#' @param model A `SSModelLeadingIndicator` object
#' 
#' @export
plot.SSModelLeadingIndicator <- function(model,title=NULL, series.name.lead="Leading Indicator", 
                                         series.name.target="Target Variable",take.log=TRUE) {
  # Call the object's plot() method if it exists
  if (!is.null(model$plot) && is.function(model$plot)) {
    model$plot(title, series.name.lead, 
               series.name.target,take.log)
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
#' 
#' @export
plot.SSModelDynamicGompertz <- function(model,title=NULL, 
                                        series.name="target variable", MA=TRUE) {
  # Call the object's plot() method if it exists
  if (!is.null(model$plot) && is.function(model$plot)) {
    model$plot(title, series.name, MA)
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

