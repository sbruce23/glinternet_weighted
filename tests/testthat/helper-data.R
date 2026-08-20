make_data <- function(n=48) {
  set.seed(1401)
  X <- cbind(sample(0:2, n, TRUE), sample(0:1, n, TRUE), rnorm(n), rnorm(n))
  y <- 0.7 + X[, 3] - 0.8*X[, 4] + 1.2*X[, 3]*X[, 4] + rnorm(n, sd=.4)
  list(X=X, y=y, levels=c(3,2,1,1), w=runif(n, .1, 2))
}
