test_that("tight Gaussian fits satisfy active-group vector KKT conditions", {
  set.seed(20260825)
  n <- 120
  X <- cbind(rnorm(n), rnorm(n), rnorm(n))
  y <- 0.4 + X[, 1] - 0.7 * X[, 2] +
    1.2 * X[, 1] * X[, 2] + rnorm(n, sd = 0.35)
  weights <- glinternet:::validate_weights(runif(n, 0.4, 2), n)
  Z <- apply(X, 2, glinternet:::standardize, weights = weights)

  active <- list(
    cat = NULL,
    cont = NULL,
    catcat = NULL,
    contcont = matrix(c(1, 2), nrow = 1),
    catcont = NULL
  )
  lambda <- 0.015
  intercept <- glinternet:::initial_intercept(y, weights, "gaussian")
  beta <- glinternet:::initialize_betahat(active, NULL, intercept, NULL)
  solution <- glinternet:::group_lasso(
    X = NULL,
    Z = Z,
    Y = y,
    weights = weights,
    activeSet = active,
    betahat = beta,
    numLevels = NULL,
    lambda = lambda,
    family = "gaussian",
    tol = 1e-9,
    maxIter = 20000,
    verbose = FALSE
  )

  product <- Z[, 1] * Z[, 2]
  product <- product - sum(weights * product) / n
  product_norm <- sqrt(sum(weights * product^2))
  interaction_design <- cbind(
    Z[, 1] / sqrt(3),
    Z[, 2] / sqrt(3),
    product / (sqrt(3) * product_norm)
  )
  gradient <- -drop(crossprod(
    interaction_design,
    weights * solution$res
  )) / n
  native_beta <- solution$betahat[-1]
  vector_kkt_residual <- sqrt(sum((gradient +
    lambda * native_beta / sqrt(sum(native_beta^2)))^2))

  expect_true(solution$converged)
  expect_lte(abs(sum(weights * solution$res) / n), 2e-8)
  expect_equal(solution$activeSet$contcont, matrix(c(1, 2), nrow = 1))
  expect_lte(vector_kkt_residual, 2e-6)
})

test_that("Gaussian integer replication uses the row-geometry lambda adjustment", {
  set.seed(20260826)
  n <- 45
  X <- cbind(rnorm(n), rnorm(n))
  y <- 0.2 + X[, 1] - 0.4 * X[, 2] +
    0.9 * X[, 1] * X[, 2] + rnorm(n, sd = 0.2)
  weights <- sample(1:4, n, replace = TRUE)
  lambda <- 0.025

  weighted <- glinternet(
    X, y, c(1, 1), lambda = lambda, weights = weights,
    tol = 1e-9, maxIter = 12000
  )
  index <- rep(seq_len(n), weights)
  N <- length(index)
  replicated <- glinternet(
    X[index, , drop = FALSE], y[index], c(1, 1),
    lambda = lambda * sqrt(n / N), tol = 1e-9, maxIter = 12000
  )

  expect_equal(
    predict(weighted, X, "response"),
    predict(replicated, X, "response"),
    tolerance = 3e-6
  )
})
