standardize = function(x, weights=rep(1, length(x))){
  positiveWeight = weights > 0
  result = x - weighted_mean(x[positiveWeight], weights[positiveWeight])
  scale = sqrt(sum(weights[positiveWeight] * result[positiveWeight]^2))
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) return(rep(0, length(x)))
  result / scale
}
