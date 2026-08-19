validate_weights = function(weights, n) {
  if (is.null(weights)) weights = rep(1, n)
  if (!is.numeric(weights) || length(weights) != n || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("weights must be a finite, nonnegative numeric vector of length n with positive sum")
  }
  total = sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop("weights must be a finite, nonnegative numeric vector of length n with positive sum")
  }
  as.numeric(weights) / total * n
}

weighted_mean = function(x, weights) {
  positiveWeight = weights > 0
  scale = max(abs(x[positiveWeight]))
  if (scale == 0) return(0)
  scaledMean = sum(weights[positiveWeight] * (x[positiveWeight] / scale)) / sum(weights[positiveWeight])
  scale * max(-1, min(1, scaledMean))
}

initial_intercept = function(Y, weights, family) {
  mu = weighted_mean(Y, weights)
  if (family == "gaussian") return(mu)
  positiveMass = sum(weights[Y == 1])
  negativeMass = sum(weights[Y == 0])
  if (!is.finite(positiveMass) || !is.finite(negativeMass) ||
      positiveMass <= 0 || negativeMass <= 0) {
    stop("weighted binomial response must contain positive weight in both outcome classes")
  }
  log(positiveMass) - log(negativeMass)
}
