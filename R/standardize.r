standardize = function(x, weights=rep(1, length(x))){
  result = x - weighted_mean(x, weights)
  scale = sqrt(sum(weights * result^2))
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) return(rep(0, length(x)))
  result / scale
}
