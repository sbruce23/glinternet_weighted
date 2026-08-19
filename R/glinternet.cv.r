glinternet.cv = function(X, Y, numLevels, nFolds=10, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, interactionCandidates=NULL, interactionPairs=NULL, screenLimit=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000, verbose=FALSE, numCores=1, weights=NULL, foldid=NULL) {

  # get call and family
  thisCall = match.call()
  family = match.arg(family)

  # make sure inputs are valid
  n = length(Y)
  rawWeights = if (is.null(weights)) rep(1, n) else weights
  validate_weights(rawWeights, n)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  stopifnot(n==nrow(X), pCat+pCont==ncol(X), family=="gaussian"||family=="binomial")

  fullfitted = glinternet(X=X, Y=Y, numLevels=numLevels, lambda=lambda, nLambda=nLambda,
                         lambdaMinRatio=lambdaMinRatio, interactionCandidates=interactionCandidates,
                         interactionPairs=interactionPairs, screenLimit=screenLimit, family=family,
                         tol=tol, maxIter=maxIter, verbose=verbose, numCores=numCores, weights=rawWeights)
  if(verbose) {
    cat("\n Done fit on all data\n")
  }

  lambda=fullfitted$lambda
  nlambda=length(lambda)

  # create the folds
  if (is.null(foldid)) {
    if (length(nFolds) != 1 || !is.finite(nFolds) || nFolds != as.integer(nFolds) || nFolds < 2 || nFolds > n) {
      stop("nFolds must be an integer between 2 and n")
    }
    folds = sample(rep(seq_len(nFolds), length.out=n))
  } else {
    if (length(foldid) != n || any(!is.finite(foldid)) || any(foldid != as.integer(foldid))) {
      stop("foldid must contain one finite integer label per observation")
    }
    folds = match(foldid, unique(foldid))
    nFolds = max(folds)
    if (nFolds < 2) stop("foldid must define at least two nonempty folds")
  }
  foldMass = vapply(seq_len(nFolds), function(k) sum(rawWeights[folds == k]), numeric(1))
  if (any(foldMass <= 0)) stop("each validation fold must have positive total weight")

  # helper for loss calculation
  compute_loss = function(y, yhat, weights, family) {
    if (family == "gaussian") {
      return (sum(weights*(y-yhat)^2)/(2*sum(weights)))
    }
    yhat = sapply(yhat, function(x) min(max(1e-15, x), 1-1e-15))
    -sum(weights * (y*log(yhat) + (1-y)*log(1-yhat)))/sum(weights)
  }
  loss = matrix(0, nFolds, nlambda)

  X=as.matrix(X)
  for (fold in 1:nFolds) {
    testIndex= (folds == fold)
    trainIndex = !testIndex
    fitted = glinternet(X=X[trainIndex,,drop=FALSE], Y=Y[trainIndex], numLevels=numLevels,
                        lambda=lambda, nLambda=nlambda, lambdaMinRatio=lambdaMinRatio,
                        interactionCandidates=interactionCandidates, interactionPairs=interactionPairs,
                        screenLimit=screenLimit, numToFind=NULL, family=family, tol=tol,
                        maxIter=maxIter, verbose=verbose, numCores=numCores,
                        weights=rawWeights[trainIndex])
    YtestHat = predict(fitted, X[testIndex,,drop=FALSE], "response")
    loss[fold, ] = apply(YtestHat, 2, function(yhat) compute_loss(Y[testIndex], yhat, rawWeights[testIndex], family))
    if(verbose) {
      cat("\n Done fold",fold,"\n")
    }
  }

  # compute cv errors and get minimum
  foldShare = foldMass / sum(foldMass)
  cv = drop(crossprod(foldShare, loss))
  effectiveFolds = 1 / sum(foldShare^2)
  foldVariance = colSums(foldShare * sweep(loss, 2, cv)^2) / (1-sum(foldShare^2))
  cvStd = sqrt(foldVariance / effectiveFolds)
  bestIndex1Std = which(cv <= min(cv)+cvStd[which.min(cv)])
  bestIndex = which.min(cv)
  lambdaHat1Std = lambda[bestIndex1Std[1]]
  lambdaHat = lambda[bestIndex]

  # return fit on full dataset with chosen lambda
  output = list(call=thisCall, glinternetFit=fullfitted, fitted=fullfitted$fitted[, bestIndex], activeSet=fullfitted$activeSet[bestIndex], betahat=fullfitted$betahat[bestIndex], lambda=lambda, lambdaHat=lambdaHat, lambdaHat1Std=lambdaHat1Std, cvErr=cv, cvErrStd=cvStd, family=family, numLevels=numLevels, nFolds=nFolds, foldid=folds, foldLoss=loss, foldWeight=foldMass, weights=rawWeights)
  class(output) = "glinternet.cv"

  return (output)
}
