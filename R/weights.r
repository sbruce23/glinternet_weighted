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

weighted_mean = function(x, weights) sum(weights * x) / sum(weights)

initial_intercept = function(Y, weights, family) {
  mu = weighted_mean(Y, weights)
  if (family == "gaussian") return(mu)
  if (mu <= 0 || mu >= 1) {
    stop("weighted binomial response must contain positive weight in both outcome classes")
  }
  qlogis(mu)
}
