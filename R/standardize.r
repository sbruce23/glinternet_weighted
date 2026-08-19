standardize = function(x, weights=rep(1, length(x))){
  positiveWeight = weights > 0
  magnitude = max(abs(x[positiveWeight]))
  if (magnitude == 0) return(rep(0, length(x)))
  scaled = x / magnitude
  center = sum(weights[positiveWeight] * scaled[positiveWeight]) / sum(weights[positiveWeight])
  result = scaled - center
  scale = sqrt(sum(weights[positiveWeight] * result[positiveWeight]^2))
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) return(rep(0, length(x)))
  result / scale
}
