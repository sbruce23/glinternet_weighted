glinternet = function(X, Y, numLevels, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, interactionCandidates=NULL, interactionPairs=NULL, screenLimit=NULL, numToFind=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000, verbose=FALSE, numCores=1, weights=NULL) {

  # get call and family
  thisCall = match.call()
  family = match.arg(family)

  # make sure inputs are valid
  X = as.matrix(X)
  n = length(Y)
  if (!is.numeric(Y) || !length(Y) || any(!is.finite(Y)))
    stop("Y must be a finite numeric vector")
  if (!is.numeric(X) || !ncol(X) || any(!is.finite(X)))
    stop("X must be a finite numeric matrix with at least one column")
  if (!is.numeric(numLevels) || length(numLevels) != ncol(X) ||
      any(!is.finite(numLevels)) || any(numLevels != as.integer(numLevels)) ||
      any(numLevels < 1))
    stop("numLevels must contain one positive integer per column of X")
  if (nrow(X) != n) stop("X and Y must contain the same number of observations")
  if (length(nLambda) != 1 || !is.finite(nLambda) || nLambda != as.integer(nLambda) || nLambda < 1)
    stop("nLambda must be a positive integer")
  if (length(lambdaMinRatio) != 1 || !is.finite(lambdaMinRatio) ||
      lambdaMinRatio <= 0 || lambdaMinRatio > 1)
    stop("lambdaMinRatio must be in (0, 1]")
  if (!is.null(screenLimit) && (length(screenLimit) != 1 || !is.finite(screenLimit) ||
      screenLimit != as.integer(screenLimit) || screenLimit < 1))
    stop("screenLimit must be NULL or a positive integer")
  if (!is.null(numToFind) && (length(numToFind) != 1 || !is.finite(numToFind) ||
      numToFind != as.integer(numToFind) || numToFind < 1))
    stop("numToFind must be NULL or a positive integer")
  if (length(tol) != 1 || !is.finite(tol) || tol <= 0)
    stop("tol must be a finite positive number")
  if (length(maxIter) != 1 || !is.finite(maxIter) || maxIter != as.integer(maxIter) || maxIter < 1)
    stop("maxIter must be a positive integer")
  if (length(numCores) != 1 || !is.finite(numCores) || numCores != as.integer(numCores) || numCores < 1)
    stop("numCores must be a positive integer")
  weights = validate_weights(weights, n)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  if (family=="binomial" && !all(Y %in% 0:1)) {
    stop("Error:family=binomial but Y not in {0,1}")
  }
  for (i in seq_len(ncol(X))) {
    if (numLevels[i]>1 && any(X[,i] != as.integer(X[,i]) | X[,i] < 0 | X[,i] >= numLevels[i])) {
      stop(sprintf("Column %d of X is categorical, but not coded as {0, 1, ...}. Refer to glinternet help on what the X argument should be.", i))
    }
  }

  contIndices = which(numLevels == 1)
  catIndices = which(numLevels > 1)
  if (!is.null(interactionCandidates)) {
    if (!is.numeric(interactionCandidates) || any(!is.finite(interactionCandidates)) ||
        any(interactionCandidates != as.integer(interactionCandidates)) ||
        any(interactionCandidates < 1 | interactionCandidates > ncol(X)))
      stop("interactionCandidates must contain valid integer column indices")
    interactionCandidates = unique(as.integer(interactionCandidates))
  }

  # specific interaction pairs
  if (!is.null(interactionPairs)) {
    # sanity check
    if (!is.matrix(interactionPairs) || ncol(interactionPairs) != 2 || !nrow(interactionPairs) ||
        !is.numeric(interactionPairs) || any(!is.finite(interactionPairs)) ||
        any(interactionPairs != as.integer(interactionPairs)) ||
        any(interactionPairs < 1 | interactionPairs > ncol(X)) ||
        any(interactionPairs[,1] == interactionPairs[,2])) {
      stop("interactionPairs must be a nonempty two-column matrix of distinct valid integer column indices")
    }
    if (!is.null(interactionCandidates)) {
      stop("If interactionPairs is set, interactionCandidates must be NULL.")
    }
    pairs = list(contcont=NULL, catcat=NULL, catcont=NULL)
    for (i in 1:nrow(interactionPairs)) {
      left = interactionPairs[i, 1]
      right = interactionPairs[i, 2]
      if (numLevels[left] == 1 && numLevels[right] == 1) {
        pairs$contcont = c(pairs$contcont, which(contIndices %in% c(left, right)))
      } else if (numLevels[left] == 1) {
        pairs$catcont = c(pairs$catcont, which(catIndices == right), which(contIndices == left))
      } else if (numLevels[right] == 1) {
        pairs$catcont = c(pairs$catcont, which(catIndices == left), which(contIndices == right))
      } else {
        pairs$catcat = c(pairs$catcat, which(catIndices %in% c(left, right)))
      }
    }
    # convert to matrices
    pairs = lapply(pairs, function(x) {
      if (!is.null(x)) {
        return(matrix(x, ncol=2, byrow=TRUE))
      } else {
        return(NULL)
      }
    })
    interactionPairs = pairs
  }

  # separate into categorical and continuous parts
  if (pCont > 0) {
    continuousCandidates = NULL
    Z = as.matrix(apply(as.matrix(X[, contIndices, drop=FALSE]), 2, standardize, weights=weights))
    if (!is.null(interactionCandidates)) {
      continuousCandidates = which(contIndices %in% interactionCandidates)
    }
  } else {
    Z = NULL
    continuousCandidates = NULL
  }
  if (pCat > 0){
    categoricalCandidates = NULL
    levels = numLevels[catIndices]
    Xcat = as.matrix(X[, catIndices])
    if (!is.null(interactionCandidates)) {
      categoricalCandidates = which(catIndices %in% interactionCandidates)
    }
  } else {
    levels = NULL
    Xcat = NULL
    categoricalCandidates = NULL
  }

  # compute variable norms
  intercept = initial_intercept(Y, weights, family)
  responseMean = if (family == "gaussian") intercept else plogis(intercept)
  res = Y - responseMean
  candidates = get_candidates(Xcat, Z, res, weights, n, pCat, pCont, levels, interactionPairs, categoricalCandidates, continuousCandidates, screenLimit, numCores=numCores)

  # lambda grid if not user provided
  userLambda = !is.null(lambda)
  if (!userLambda) {
    lambda = get_lambda_grid(candidates, nLambda, lambdaMinRatio)
    lambdaMax = lambda[1]
  } else {
    if (!is.numeric(lambda) || !length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
      stop("lambda must be a finite positive numeric vector")
    }
    if (any(diff(lambda) > 0)) {
      stop("Error: input lambda sequence is not monotone decreasing.")
    }
    lambdaMax = max(vapply(candidates$norms, function(x)
      if (is.null(x) || !length(x)) 0 else max(x), numeric(1)))
    nLambda = length(lambda)
  }

  # initialize storage for results
  fitted = matrix(responseMean, n, nLambda)
  activeSet = vector("list", nLambda)
  betahat = vector("list", nLambda)
  objValue = rep(0, nLambda)
  converged = rep(TRUE, nLambda)
  iterations = integer(nLambda)
  softplus = pmax(intercept, 0) + log1p(exp(-abs(intercept)))
  positiveWeight = weights > 0
  nullObjective = ifelse(family=="gaussian",
                         sum(weights[positiveWeight] * res[positiveWeight]^2)/(2*n),
                         sum(weights[positiveWeight] * (softplus - Y[positiveWeight] * intercept))/n)
  if (!userLambda) {
    betahat[[1]] = intercept
    objValue[1] = nullObjective
  }
  degeneratePath = !userLambda && lambdaMax <= 0
  if (degeneratePath) {
    emptyActive = setNames(vector("list",5),c("cat","cont","catcat","contcont","catcont"))
    for (j in seq_len(nLambda)) {
      activeSet[[j]] = emptyActive
      betahat[[j]] = intercept
      objValue[j] = nullObjective
    }
  }

  # ever-active set + sequential strong rules + group lasso
  loopStart = if (degeneratePath) nLambda + 1 else if (userLambda) 1 else 2
  fitIndices = if (loopStart <= nLambda) seq.int(loopStart, nLambda) else integer()
  lastIndex = if (degeneratePath) nLambda else if (userLambda) 0 else 1
  for (i in fitIndices){
    lastIndex = i
    if (verbose) {
      cat("lambda ", i, ": ", lambda[i], "\n")
    }
    if (i == 1) {
      previousActive = setNames(vector("list", 5), c("cat", "cont", "catcat", "contcont", "catcont"))
      previousBeta = intercept
      previousLambda = lambdaMax
    } else {
      previousActive = activeSet[[i-1]]
      previousBeta = betahat[[i-1]]
      previousLambda = lambda[i-1]
    }
    activeSet[[i]] = strong_rules(candidates, lambda[i], previousLambda)
    betahat[[i]] = initialize_betahat(activeSet[[i]], previousActive, previousBeta, levels)
    while (TRUE) {
      # group lasso on strong set
      solution = group_lasso(Xcat, Z, Y, weights, activeSet[[i]], betahat[[i]], levels, lambda[i], family, tol, maxIter, verbose)
      activeSet[[i]] = solution$activeSet
      betahat[[i]] = solution$betahat
      res = solution$res
      objValue[i] = solution$objValue
      converged[i] = solution$converged
      iterations[i] = iterations[i] + solution$iterations
      if (!solution$converged) {
        warning(sprintf("FISTA reached maxIter=%d without convergence at lambda index %d", maxIter, i),
                call.=FALSE)
        break
      }
      # check kkt conditions on the rest
      check = check_kkt(Xcat, Z, res, weights, n, pCat, pCont, levels, candidates, activeSet[[i]], lambda[i], numCores)
      candidates$norms = check$norms
      if (check$flag) {
        break
      }
      betahat[[i]] = initialize_betahat(check$activeSet, activeSet[[i]], betahat[[i]], levels)
      activeSet[[i]] = check$activeSet
    }
    # update the candidate set if necessary
    if (!is.null(screenLimit) && (screenLimit<pCat+pCont) && i<nLambda) {
      candidates = get_candidates(Xcat, Z, res, weights, n, pCat, pCont, levels, interactionPairs, categoricalCandidates, continuousCandidates, screenLimit, activeSet[[i]], candidates$norms, numCores)
    }
    # get fitted values
    fitted[, i] = Y - res
    # compute total number of interactions found
    if (!is.null(numToFind)) {
      numFound = sum(sapply(activeSet[[i]][3:5], function(x) ifelse(is.null(x), 0, nrow(x))))
      if (numFound >= numToFind) {
        break
      }
    }
  }

  # rescale betahat
  i = lastIndex
  Z = as.matrix(X[, numLevels==1])
  betahatRescaled = lapply(1:i, function(j) rescale_betahat(activeSet[[j]], betahat[[j]], Xcat, Z, weights, levels, n))

  output = list(call=thisCall, fitted=fitted[, 1:i, drop=FALSE], lambda=lambda[1:i], objValue=objValue[1:i],
                activeSet=activeSet[1:i], betahat=betahatRescaled[1:i], numLevels=numLevels,
                family=family, weights=weights, converged=converged[1:i],
                iterations=iterations[1:i])
  class(output) = "glinternet"

  return (output)
}
