test_that("weights are validated and normalized", {
  validate <- glinternet:::validate_weights
  expect_equal(validate(NULL, 3), rep(1, 3))
  expect_equal(sum(validate(c(1, 2, 3), 3)), 3)
  expect_equal(validate(c(1, 2, 3), 3), validate(c(10, 20, 30), 3))
  expect_error(validate(c(1, 2), 3), "length n")
  expect_error(validate(c(1, -1, 2), 3), "nonnegative")
  expect_error(validate(c(1, Inf, 2), 3), "finite")
  expect_error(validate(c(0, 0, 0), 3), "positive sum")
  expect_error(validate(rep(.Machine$double.xmax, 3), 3), "finite")
  expect_error(validate(c("1", "2", "3"), 3), "numeric vector")
  expect_true(all(is.finite(validate(c(.Machine$double.xmax/2, 1, 1), 3))))
})

test_that("weighted standardization has exact weighted geometry", {
  x <- c(-2, 0, 4, 1000)
  v <- glinternet:::validate_weights(c(1, 2, 3, 0), length(x))
  z <- glinternet:::standardize(x, v)
  expect_equal(sum(v * z), 0, tolerance=1e-14)
  expect_equal(sum(v * z^2), 1, tolerance=1e-14)
  expect_equal(z[1:3], glinternet:::standardize(x[1:3], v[1:3]), tolerance=1e-14)
})

test_that("all-one weights preserve legacy results and weight scale is irrelevant", {
  d <- make_data()
  args <- list(X=d$X, Y=d$y, numLevels=d$levels, nLambda=5, maxIter=3000)
  f0 <- do.call(glinternet, args)
  f1 <- do.call(glinternet, c(args, list(weights=rep(1, length(d$y)))))
  fw <- do.call(glinternet, c(args, list(weights=d$w)))
  fws <- do.call(glinternet, c(args, list(weights=73*d$w)))
  expect_equal(f0$lambda, f1$lambda, tolerance=1e-14)
  expect_equal(f0$fitted, f1$fitted, tolerance=1e-12)
  expect_equal(fw$lambda, fws$lambda, tolerance=1e-12)
  expect_equal(fw$fitted, fws$fitted, tolerance=1e-6)
  expect_equal(fw$objValue, fws$objValue, tolerance=1e-6)
})

test_that("weighted fits are finite and original-scale predictions agree", {
  d <- make_data()
  fg <- glinternet(d$X, d$y, d$levels, nLambda=5, weights=d$w, maxIter=3000)
  yb <- as.numeric(d$y > median(d$y))
  fb <- glinternet(d$X, yb, d$levels, nLambda=4, weights=d$w,
                   family="binomial", maxIter=3000)
  expect_true(all(is.finite(fg$objValue)))
  expect_true(all(is.finite(fb$objValue)))
  expect_equal(predict(fg, d$X, "response"), fg$fitted, tolerance=1e-7)
  expect_equal(predict(fb, d$X, "response"), fb$fitted, tolerance=1e-7)
})

test_that("weighted CV aggregates losses and uncertainty by fold mass", {
  d <- make_data(36)
  foldid <- rep(1:3, each=12)
  fit <- glinternet.cv(d$X, d$y, d$levels, nFolds=3, nLambda=4,
                       weights=d$w, foldid=foldid, maxIter=2000)
  a <- fit$foldWeight/sum(fit$foldWeight)
  direct <- drop(crossprod(a, fit$foldLoss))
  keff <- 1/sum(a^2)
  variance <- colSums(a*sweep(fit$foldLoss, 2, direct)^2)/(1-sum(a^2))
  expect_equal(fit$cvErr, direct, tolerance=1e-12)
  expect_equal(fit$cvErrStd, sqrt(variance/keff), tolerance=1e-12)
  eligible <- which(fit$cvErr <= min(fit$cvErr)+fit$cvErrStd[which.min(fit$cvErr)])
  expect_equal(fit$lambdaHat1Std, fit$lambda[eligible[1]])
})

test_that("cross-validation rejects invalid fold definitions", {
  d <- make_data(12)
  expect_error(glinternet.cv(d$X,d$y,d$levels,nFolds=1,weights=d$w),"between 2 and n")
  expect_error(glinternet.cv(d$X,d$y,d$levels,nFolds=13,weights=d$w),"between 2 and n")
  expect_error(glinternet.cv(d$X,d$y,d$levels,weights=d$w,foldid=rep(1,12)),"at least two")
  expect_error(glinternet.cv(d$X,d$y,d$levels,weights=c(rep(1,6),rep(0,6)),
                             foldid=rep(1:2,each=6)),"positive total weight")
  expect_error(glinternet.cv(d$X,d$y,d$levels,weights=d$w,foldid=c(rep(1,11),1.5)),
               "integer label")
})

test_that("zero-weight observations do not affect weighted moments", {
  x1 <- c(1,2,3,999)
  x2 <- c(1,2,3,-999)
  v <- glinternet:::validate_weights(c(1,2,1,0), 4)
  expect_equal(glinternet:::standardize(x1, v)[1:3],
               glinternet:::standardize(x2, v)[1:3], tolerance=1e-14)
})

test_that("degenerate continuous products stay finite", {
  set.seed(44)
  n <- 40
  x <- rep(c(-1, 1), each=n/2)
  X <- cbind(x, x, rnorm(n))
  y <- 1 + x + rnorm(n, sd=.2)
  fit <- glinternet(X, y, rep(1, 3), nLambda=5, weights=seq_len(n), maxIter=3000)
  expect_true(all(is.finite(fit$objValue)))
  expect_true(all(is.finite(fit$fitted)))
  expect_equal(predict(fit, X, "response"), fit$fitted, tolerance=1e-7)
})
