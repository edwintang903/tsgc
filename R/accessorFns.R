# Trimmed accessor functions for the idx_series refactor smoke test.
# These are the plain (non S3-dispatched) accessor helpers from the
# original accessorFns.R that have no dependency on xts/dates and are
# still valid as-is. The S3 print/summary/plot wrappers for
# SSModelDynamicGompertz/FilterResults/etc. are NOT included here because
# plotting has been removed in this pass of the refactor; print/summary are
# called directly as RefClass methods in the test script instead.

#' @export
output<-function(object){
  return(object$output)
}

#' @export
modelKFS<-function(object){
  return(object$model)
}

#' @export
seasonalComp<-function(object){
  attr(modelKFS(object)$terms, "specials")$SSMseasonal
}

#' @export
att<-function(object){
  object$att
}

#' @export
Ptt<-function(object){
  object$Ptt
}

#' @export
get_V<-function(object){
  object$V
}

#' @export
matrixKFS<-function(object,matrix){
  modelKFS(object)[[matrix]]
}

#' @export
gety<-function(object){
  object$y
}

#' @export
gety.hat<-function(object){
  object$y.hat
}

#' @export
alphahat<-function(object){
  object$alphahat
}

#' @export
estimate<-function(model){
  if (!inherits(model, "SSModelDynamicGompertz") && !inherits(model, "SSModelLeadingIndicator")){
    stop("model must be a SSModelDynamicGompertz or SSModelLeadingIndicator object.")
  }
  model$estimate()
}