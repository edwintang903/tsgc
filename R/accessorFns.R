#' @title Extract output of FilterResults or FilterResultsLI
#'
#' @description Accessor method to access the fitted KFS model from
#' `FilterResults` or `FilterResultsLI`.
#'
#' @param object A `FilterResults` or `FilterResultsLI` object.
#'
#' @returns The fitted KFS model.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' output(res)
#'
#' @export
output <- function(object) {
  return(object$output)
}

#' @title Extract SSModel object within a KFS object
#'
#' @description Accessor method to access the fitted `SSModel`.
#'
#' @param object A `KFS` object.
#'
#' @returns The `SSModel` object underlying the `KFS` fit.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' modelKFS(output(res))
#'
#' @export
modelKFS <- function(object) {
  return(object$model)
}

#' @title Extract number of seasonal components used in a KFS object
#'
#' @description Accessor method to access the number of seasonal components
#' used in a `KFS` object.
#'
#' @param object A `KFS` object.
#'
#' @returns The seasonal component specification stored on the model terms.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' seasonalComp(output(res))
#'
#' @export
seasonalComp <- function(object) {
  attr(modelKFS(object)$terms, "specials")$SSMseasonal
}

#' @title Extract filtered state estimates from a KFS object
#'
#' @description Accessor method to access the filtered state estimates in a
#' `KFS` object.
#'
#' @param object A `KFS` object.
#'
#' @returns The filtered state estimates.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' att(output(res))
#'
#' @export
att <- function(object) {
  object$att
}

#' @title Extract error covariance matrix of filtered states from a KFS object
#'
#' @description Accessor method to access the non-diffuse part of the error
#' covariance matrix of the filtered states in a `KFS` object.
#'
#' @param object A `KFS` object.
#'
#' @returns The error covariance matrix of the filtered states.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' Ptt(output(res))
#'
#' @export
Ptt <- function(object) {
  object$Ptt
}

#' @title Extract error covariance matrix of smoothed states from a KFS object
#'
#' @description Accessor method to access the non-diffuse part of the error
#' covariance matrix of the smoothed states in a `KFS` object.
#'
#' @param object A `KFS` object.
#'
#' @returns The error covariance matrix of the smoothed states.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' get_V(output(res))
#'
#' @export
get_V <- function(object) {
  object$V
}

#' @title Extract matrices of the observation, state, and disturbance
#' equations from a KFS object
#'
#' @description Accessor method to access the matrices used in the
#' observation, state, and disturbance equations of a `KFS` object.
#'
#' @param object A `KFS` object.
#' @param matrix Character string naming a matrix component of `SSModel`,
#' e.g. `"H"`, `"T"`, `"R"`, `"Q"`.
#'
#' @returns The requested matrix.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' matrixKFS(output(res), "Z")
#'
#' @export
matrixKFS <- function(object, matrix) {
  modelKFS(object)[[matrix]]
}

#' @title Extract time series y from an SSModel object
#'
#' @description Accessor method to access the time series `y` in an
#' `SSModel` object.
#'
#' @param object An `SSModel` object.
#'
#' @returns The time series `y`.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' gety(modelKFS(output(res)))
#'
#' @export
gety <- function(object) {
  object$y
}

#' @title Extract prediction y.hat from predict.all output
#'
#' @description Accessor method to access the prediction `y.hat` in the
#' output of `predict_all`.
#'
#' @param object An object returned by `predict_all`.
#'
#' @returns The predicted values `y.hat`.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' all_predictions <- res$predict_all(n.ahead = 7)
#' gety.hat(all_predictions)
#'
#' @export
gety.hat <- function(object) {
  object$y.hat
}

#' @title Extract alphahat from a KFS object
#'
#' @description Accessor method to access `alphahat`, the smoothed state
#' estimates, from a fitted `KFS` object.
#'
#' @param object A `KFS` object.
#'
#' @returns The smoothed state estimates `alphahat`.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' res <- estimate(model)
#' alphahat(output(res))
#'
#' @export
alphahat <- function(object) {
  object$alphahat
}

#' @title Estimate an SSModelDynamicGompertz or SSModelLeadingIndicator model
#'
#' @description Accessor method that calls the `estimate` method of a
#' `SSModelDynamicGompertz` or `SSModelLeadingIndicator` object.
#'
#' @param model A `SSModelDynamicGompertz` or `SSModelLeadingIndicator`
#' object.
#'
#' @returns A `FilterResults` or `FilterResultsLI` object containing the
#' estimation results.
#'
#' @examples
#' library(tsgc)
#' data(gauteng, package = "tsgc")
#' conv <- xts_to_idx(gauteng)
#' model <- SSModelDynamicGompertz$new(Y = conv$series, q = 0.005,
#'                                     calendar = conv$calendar,
#'                                     end = idx_to_pos(conv$calendar, "2020-07-20"))
#' estimate(model)
#'
#' @export
estimate <- function(model) {
  if (!inherits(model, "SSModelDynamicGompertz") && !inherits(model, "SSModelLeadingIndicator")) {
    stop("model must be a SSModelDynamicGompertz or SSModelLeadingIndicator object.")
  }
  model$estimate()
}