get_lambda_grid = function(candidates, nLambda, lambdaMinRatio){
  
  lambdaMax = max(sapply(candidates$norms, function(x) ifelse(is.null(x), 0, max(x))))
  if (!is.finite(lambdaMax) || lambdaMax <= 0) return(rep(0, nLambda))
  if (nLambda == 1) return(lambdaMax)
  lambdaMin = lambdaMinRatio * lambdaMax
  f = seq(0,1,1/(nLambda-1))
  lambda = lambdaMax^(1-f) * lambdaMin^f
}

