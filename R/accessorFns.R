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

#' @title Extract output of FilterResults
#
#' @description Accessor method to access the fitted KFS model
#'
#' @param object FilterResults object
#'
#'
#' @export
output<-function(object){
  return(object$output)
}

#' @title Extract SSModel of KFS
#
#' @description Accessor method to access the fitted SSModel
#'
#' @param object FilterResults object
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
#' @description Accessor method to access time series y in SSModel object
#'
#' @param object SSModel object
#'
#' @export
gety<-function(object){
  object$y
}

#' @title Extract prediction y.hat in SSModel
#
#' @description Accessor method to access time series y in SSModel object
#'
#' @param object SSModel object
#'
#' @export
gety.hat<-function(object){
  object$y.hat
}


#' @title Extract alphahat in SSModel
#
#' @description Accessor method to access alphahat in SSModel object
#'
#' @param object SSModel object
#'
#' @export
alphahat<-function(object){
  object$alphahat
}
#' 
#' #' @title Calling estimate method for SSModelDynamicGompertz class
#' #
#' #' @description Accessor method to obtain estimated model for SSModelDynamicGompertz class
#' #'
#' #' @param model SSModelDynamicGompertz object
#' #'
#' #' @export
#' estimate<-function(model){
#'   model$estimate()
#' }

