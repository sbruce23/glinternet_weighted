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
  variance <- colSums(a*sweep(fit$foldLoss, 2, direct)^2)/sum(a*(1-a))
  expect_equal(fit$cvErr, direct, tolerance=1e-12)
  expect_equal(fit$cvErrStd, sqrt(variance/keff), tolerance=1e-12)
  eligible <- which(fit$cvErr <= min(fit$cvErr)+fit$cvErrStd[which.min(fit$cvErr)])
  expect_equal(fit$lambdaHat1Std, fit$lambda[eligible[1]])
})

test_that("CV uncertainty stays finite for extremely unequal fold masses", {
  set.seed(323)
  X <- matrix(rnorm(12), 6, 2)
  y <- X[,1] + rnorm(6, sd=.1)
  folds <- rep(1:3, each=2)
  w <- c(5, 5, rep(5e-200, 4))
  fit <- glinternet.cv(X, y, c(1,1), lambda=.2, weights=w, foldid=folds,
                       maxIter=4000)
  a <- fit$foldWeight/sum(fit$foldWeight)
  direct <- drop(crossprod(a, fit$foldLoss))
  variance <- colSums(a*sweep(fit$foldLoss, 2, direct)^2)/sum(a*(1-a))
  expected <- sqrt(variance/(1/sum(a^2)))
  expect_true(all(is.finite(fit$cvErrStd)))
  expect_equal(fit$cvErrStd, expected, tolerance=1e-12)
  expect_true(is.finite(fit$lambdaHat1Std))
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

test_that("malformed model inputs fail clearly", {
  d <- make_data(12)
  expect_error(glinternet(d$X,d$y,d$levels,interactionPairs=matrix(c(0,2),1,2)),
               "distinct valid integer")
  expect_error(glinternet(d$X,d$y,d$levels,interactionPairs=matrix(c(1,1),1,2)),
               "distinct valid integer")
  expect_error(glinternet(d$X,d$y,d$levels,interactionPairs=matrix(1:3,1,3)),
               "two-column")
  expect_error(glinternet(d$X,d$y,d$levels,interactionCandidates=99),"valid integer")
  bad <- d$X; bad[1,1] <- -1
  expect_error(glinternet(bad,d$y,d$levels),"categorical")
  bad <- d$X; bad[1,1] <- .5
  expect_error(glinternet(bad,d$y,d$levels),"categorical")
  expect_error(glinternet(d$X,d$y,d$levels,lambda=c(.1,NA)),"finite positive")
  expect_error(glinternet(d$X,d$y,d$levels,lambda=c(.1,.2)),"strictly decreasing")
  expect_error(glinternet(d$X,d$y,d$levels,lambda=c(.1,.1,.05)),"strictly decreasing")
  expect_error(glinternet(d$X,d$y,d$levels,nLambda=0),"positive integer")
  expect_error(glinternet(d$X,d$y,d$levels,lambdaMinRatio=2),"lambdaMinRatio")
  expect_error(glinternet.cv(d$X,d$y,d$levels,nFolds=3,nLambda=0),"positive integer")
  expect_error(glinternet.cv(d$X,d$y,d$levels,nFolds=3,lambdaMinRatio=0),"lambdaMinRatio")
  expect_error(glinternet(d$X,d$y,d$levels,screenLimit=0),"positive integer")
  expect_error(glinternet(d$X,d$y,d$levels,numToFind=0),"positive integer")
  expect_error(glinternet(d$X,d$y,d$levels,tol=0),"finite positive")
  expect_error(glinternet(d$X,d$y,d$levels,maxIter=1.5),"positive integer")
  expect_error(glinternet(d$X,d$y,d$levels,numCores=0),"positive integer")
  expect_error(glinternet(d$X,d$y,c(.Machine$integer.max+1,1,1)),"positive integer")
  expect_error(glinternet(d$X,d$y,d$levels,nLambda=.Machine$integer.max+1),"positive integer")
  expect_error(glinternet.cv(d$X,d$y,d$levels,foldid=rep(.Machine$integer.max+1,12)),
               "integer label")
})

test_that("near-equal lambdas remain distinct and exactly selectable", {
  d <- make_data(24)
  lambda <- c(.1,.1-1e-14,.05)
  fit <- glinternet(d$X,d$y,d$levels,lambda=lambda,weights=d$w,maxIter=4000)
  expect_equal(fit$lambda,lambda,tolerance=0)
  expect_equal(predict(fit,d$X,lambda=lambda[2]),fit$fitted[,2,drop=FALSE],tolerance=1e-10)
  expect_equal(coef(fit,2)[[1]],glinternet:::extract_effects(
    fit$betahat[[2]],fit$activeSet[[2]],fit$numLevels))
})

test_that("CV folds are reproducible and normalized from arbitrary labels", {
  d <- make_data(18)
  labels <- rep(c(30,10,90),each=6)
  fit <- glinternet.cv(d$X,d$y,d$levels,lambda=c(.1,.05),weights=d$w,
                       foldid=labels,maxIter=4000)
  expect_equal(fit$foldid,rep(1:3,each=6))
  expect_equal(dim(fit$foldLoss),c(3L,2L))
  expect_equal(length(fit$foldWeight),3L)
  expect_equal(length(fit$cvErr),2L)
  expect_equal(length(fit$cvErrStd),2L)
  expect_equal(length(fit$fitted),18L)
  expect_equal(dim(predict(fit,d$X)),c(18L,1L))

  set.seed(326)
  a <- glinternet.cv(d$X,d$y,d$levels,nFolds=length(d$y),lambda=.1,
                     weights=d$w,maxIter=4000)
  set.seed(326)
  b <- glinternet.cv(d$X,d$y,d$levels,nFolds=length(d$y),lambda=.1,
                     weights=d$w,maxIter=4000)
  expect_equal(a$foldid,b$foldid)
  expect_equal(a$cvErr,b$cvErr,tolerance=1e-12)
  expect_equal(dim(a$foldLoss),c(18L,1L))
})

test_that("interaction pairs are unordered and deduplicated", {
  set.seed(325)
  X <- matrix(rnorm(120),40,3)
  y <- X[,1]+X[,2]+2*X[,1]*X[,2]+rnorm(40,sd=.1)
  canonical <- glinternet(X,y,rep(1,3),lambda=.01,
                          interactionPairs=matrix(c(1,2),1,2),maxIter=5000)
  repeated <- glinternet(X,y,rep(1,3),lambda=.01,
                         interactionPairs=rbind(c(1,2),c(2,1),c(1,2)),maxIter=5000)
  expect_equal(repeated$fitted,canonical$fitted,tolerance=1e-10)
  expect_equal(repeated$objValue,canonical$objValue,tolerance=1e-10)
  expect_equal(repeated$activeSet,canonical$activeSet)
  expect_equal(coef(repeated),coef(canonical),tolerance=1e-10)
})

test_that("length-one lambda works through weighted cross-validation", {
  d <- make_data(24)
  folds <- rep(1:3,each=8)
  fit <- glinternet.cv(d$X,d$y,d$levels,lambda=.02,weights=d$w,foldid=folds,
                       maxIter=4000)
  expect_equal(length(fit$lambda),1)
  expect_equal(dim(fit$glinternetFit$fitted),c(24L,1L))
  expect_equal(dim(predict(fit$glinternetFit,d$X)),c(24L,1L))
  expect_equal(length(fit$cvErr),1)
})

test_that("early interaction stopping keeps every path component aligned", {
  set.seed(324)
  X <- matrix(rnorm(160), 40, 4)
  y <- X[,1] + X[,2] + 3*X[,1]*X[,2] + rnorm(40, sd=.05)
  requested <- c(.2, .1, .05, .02, .01)
  fit <- glinternet(X, y, rep(1,4), lambda=requested, numToFind=1, maxIter=5000)
  pathLength <- length(fit$lambda)
  expect_true(pathLength < length(requested))
  expect_equal(ncol(fit$fitted), pathLength)
  expect_equal(length(fit$objValue), pathLength)
  expect_equal(length(fit$activeSet), pathLength)
  expect_equal(length(fit$betahat), pathLength)
  expect_equal(length(fit$converged), pathLength)
  expect_equal(length(fit$iterations), pathLength)
  expect_equal(ncol(predict(fit, X)), pathLength)
  expect_equal(length(coef(fit)), pathLength)
})

test_that("CV coefficient methods select the exact requested path index", {
  d <- make_data(30)
  folds <- rep(1:3, each=10)
  fit <- glinternet.cv(d$X, d$y, d$levels, nLambda=4, weights=d$w,
                       foldid=folds, maxIter=4000)
  for (kind in c("lambdaHat", "lambdaHat1Std")) {
    selected <- if (kind=="lambdaHat") fit$lambdaHat else fit$lambdaHat1Std
    idx <- match(selected, fit$glinternetFit$lambda)
    expect_equal(coef(fit, lambdaType=kind), coef(fit$glinternetFit, idx)[[1]])
  }
})

test_that("single-lambda S3 methods and categorical prediction validation are stable", {
  d <- make_data(30)
  fit <- glinternet(d$X, d$y, d$levels, lambda=.05, weights=d$w, maxIter=4000)
  expect_true(length(capture.output(print(fit))) > 0)
  expect_equal(dim(predict(fit,d$X)),c(30L,1L))
  expect_error(predict(fit,d$X,lambda=.05+1e-15),"not used")
  unseen <- d$X; unseen[1,1] <- d$levels[1]
  expect_error(predict(fit,unseen),"fitted integer codes")
  fractional <- d$X; fractional[1,1] <- .5
  expect_error(predict(fit,fractional),"fitted integer codes")

  yb <- as.numeric(d$y > median(d$y))
  binfit <- glinternet(d$X,yb,d$levels,lambda=.05,weights=d$w,
                       family="binomial",maxIter=4000)
  link <- predict(binfit,d$X,type="link")
  response <- predict(binfit,d$X,type="response")
  expect_equal(response,plogis(link),tolerance=1e-12)
})

test_that("single-lambda CV print plot and selections remain usable", {
  d <- make_data(30)
  folds <- rep(1:3,each=10)
  fit <- glinternet.cv(d$X,d$y,d$levels,lambda=.05,weights=d$w,
                       foldid=folds,maxIter=4000)
  expect_true(length(capture.output(print(fit))) > 0)
  pngfile <- tempfile(fileext=".png")
  grDevices::png(pngfile)
  plot(fit)
  grDevices::dev.off()
  expect_true(file.exists(pngfile))
  expect_equal(dim(predict(fit,d$X)),c(30L,1L))
  expect_equal(coef(fit),coef(fit$glinternetFit,1)[[1]])
})

test_that("degenerate zero-score data return a finite null path", {
  X <- matrix(1,20,3)
  y <- rep(0,20)
  fit <- glinternet(X,y,rep(1,3),nLambda=4,weights=seq_len(20))
  expect_equal(fit$lambda,rep(0,4))
  expect_equal(dim(fit$fitted),c(20L,4L))
  expect_true(all(is.finite(c(fit$objValue,fit$fitted))))
  expect_true(all(fit$converged))
})

test_that("automatic one-lambda path returns the null model with stable shape", {
  set.seed(322)
  X <- matrix(rnorm(60),20,3)
  y <- X[,1]+rnorm(20)
  fit <- glinternet(X,y,rep(1,3),nLambda=1,weights=seq_len(20))
  expect_equal(length(fit$lambda),1)
  expect_equal(dim(fit$fitted),c(20L,1L))
  expect_true(fit$converged[1])
  expect_equal(predict(fit,X),fit$fitted,tolerance=1e-12)
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

test_that("binomial null intercept resolves extreme minority class mass", {
  y <- c(1,0)
  for (exponent in c(-16,-50,-200)) {
    w <- glinternet:::validate_weights(c(10^exponent,1),2)
    intercept <- glinternet:::initial_intercept(y,w,"binomial")
    expect_true(is.finite(intercept))
    expect_equal(intercept,log(w[1])-log(w[2]),tolerance=1e-14)
  }
  expect_error(glinternet:::initial_intercept(y,c(0,2),"binomial"),"both outcome classes")
})
