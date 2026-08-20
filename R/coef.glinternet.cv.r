coef.glinternet.cv = function(object, lambdaType=c("lambdaHat", "lambdaHat1Std"), ...){

  lambdaType = match.arg(lambdaType)
  selectedLambda = if (lambdaType=="lambdaHat") object$lambdaHat else object$lambdaHat1Std
  idx = match(selectedLambda, object$glinternetFit$lambda)
  if (is.na(idx)) stop("Selected lambda is not present in the fitted path")
  extract_effects(object$glinternetFit$betahat[[idx]], object$glinternetFit$activeSet[[idx]], object$numLevels)
}
