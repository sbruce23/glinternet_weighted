test_that("all-one weights reproduce legacy glinternet and CV fits", {
  set.seed(20260820)
  n <- 36
  X <- cbind(rnorm(n), rnorm(n), rnorm(n))
  y <- 0.8 + 0.7*X[,1] - 0.4*X[,2] + 0.9*X[,1]*X[,2] +
    rnorm(n, sd=0.15)
  numLevels <- rep(1, ncol(X))
  lambda <- c(0.12, 0.06, 0.03)
  foldid <- rep(1:3, length.out=n)

  # Omitting weights is the legacy glinternet interface. Supplying all-one
  # weights must reproduce that fit, including the full regularization path.
  legacyFit <- glinternet(X, y, numLevels, lambda=lambda,
                          maxIter=5000, tol=1e-7)
  equalWeightFit <- glinternet(X, y, numLevels, lambda=lambda,
                               weights=rep(1, n), maxIter=5000, tol=1e-7)

  expect_equal(equalWeightFit$lambda, legacyFit$lambda, tolerance=0)
  expect_equal(equalWeightFit$fitted, legacyFit$fitted, tolerance=1e-12)
  expect_equal(equalWeightFit$objValue, legacyFit$objValue, tolerance=1e-12)
  expect_equal(equalWeightFit$activeSet, legacyFit$activeSet)
  expect_equal(equalWeightFit$betahat, legacyFit$betahat, tolerance=1e-12)
  expect_equal(predict(equalWeightFit, X), predict(legacyFit, X),
               tolerance=1e-12)

  # Fixed folds make the CV comparison deterministic and test both the
  # validation summaries and the selected full-data fits.
  legacyCv <- glinternet.cv(X, y, numLevels, lambda=lambda, foldid=foldid,
                            maxIter=5000, tol=1e-7)
  equalWeightCv <- glinternet.cv(X, y, numLevels, lambda=lambda,
                                 weights=rep(1, n), foldid=foldid,
                                 maxIter=5000, tol=1e-7)

  expect_equal(equalWeightCv$foldLoss, legacyCv$foldLoss, tolerance=1e-12)
  expect_equal(equalWeightCv$foldWeight, legacyCv$foldWeight, tolerance=0)
  expect_equal(equalWeightCv$cvErr, legacyCv$cvErr, tolerance=1e-12)
  expect_equal(equalWeightCv$cvErrStd, legacyCv$cvErrStd, tolerance=1e-12)
  expect_equal(equalWeightCv$lambdaHat, legacyCv$lambdaHat, tolerance=0)
  expect_equal(equalWeightCv$lambdaHat1Std, legacyCv$lambdaHat1Std,
               tolerance=0)
  expect_equal(predict(equalWeightCv, X, lambdaType="lambdaHat"),
               predict(legacyCv, X, lambdaType="lambdaHat"),
               tolerance=1e-12)
  expect_equal(predict(equalWeightCv, X, lambdaType="lambdaHat1Std"),
               predict(legacyCv, X, lambdaType="lambdaHat1Std"),
               tolerance=1e-12)
})

test_that("one dominant observation strongly controls a weighted fit", {
  nOrdinary <- 20
  X <- matrix(c(seq(-1, 1, length.out=nOrdinary), 8), ncol=1)
  y <- c(seq(-0.2, 0.2, length.out=nOrdinary), 100)
  dominant <- length(y)
  weights <- c(rep(1e-6, nOrdinary), 1e6)

  # A deliberately large lambda leaves only the unpenalized intercept. This
  # isolates the effect of observation weights from variable selection.
  equalWeightFit <- glinternet(X, y, 1, lambda=100,
                               maxIter=5000, tol=1e-9)
  dominantWeightFit <- glinternet(X, y, 1, lambda=100, weights=weights,
                                  maxIter=5000, tol=1e-9)

  expectedWeightedLevel <- weighted.mean(y, weights)
  equalWeightLevel <- mean(y)
  weightedError <- abs(dominantWeightFit$fitted[dominant,1] - y[dominant])
  equalWeightError <- abs(equalWeightFit$fitted[dominant,1] - y[dominant])

  expect_gt(weights[dominant]/sum(weights), 1 - 1e-10)
  expect_true(all(vapply(dominantWeightFit$activeSet[[1]], is.null,
                         logical(1))))
  expect_true(all(vapply(equalWeightFit$activeSet[[1]], is.null,
                         logical(1))))
  expect_equal(drop(dominantWeightFit$fitted),
               rep(expectedWeightedLevel, length(y)), tolerance=1e-10)
  expect_equal(drop(equalWeightFit$fitted),
               rep(equalWeightLevel, length(y)), tolerance=1e-12)
  expect_lt(weightedError, 1e-6)
  expect_gt(equalWeightError, 90)
  expect_lt(weightedError, equalWeightError*1e-8)
})
