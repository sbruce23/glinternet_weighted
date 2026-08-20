group_lasso = function(X, Z, Y, weights, activeSet, betahat, numLevels, lambda, family, tol, maxIter, verbose){

  #get size of each group, number of groups in each category, and total number of groups
  groupSizes = get_group_sizes(activeSet, numLevels)
  numGroups = sapply(activeSet, function(x) if (is.null(x)) 0 else nrow(x))
  totalGroups = sum(numGroups)

  #if active set is empty, just return estimate of the intercept
  if (totalGroups == 0){
    betahat = initial_intercept(Y, weights, family)
    mu = if (family == "gaussian") betahat else plogis(betahat)
    res = Y - mu
    positiveWeight = weights > 0
    objValue = ifelse(family=="gaussian",
                      sum(weights[positiveWeight] * res[positiveWeight]^2)/(2*length(Y)),
                      sum(weights[positiveWeight] * (pmax(betahat, 0) + log1p(exp(-abs(betahat))) - Y[positiveWeight] * betahat))/length(Y))
    return(list(betahat=betahat, activeSet=activeSet, res=res, objValue=objValue,
                converged=TRUE, iterations=0L))
  }

  n = length(Y)
  indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)

  #fit and get new betahat, res, objValue
  fit = .Call("R_gl_solver", X, Z, Y, weights, n, betahat[1], betahat[-1], numLevels, numGroups, indices$cat, indices$cont, indices$catcat, indices$contcont, indices$catcont, lambda, tol, 0.1, maxIter, ifelse(family=="gaussian", 0, 1), ifelse(verbose, 1, 0))
  res = fit$res
  objValue = fit$objValue

  #get the nonzero parts of betahat and update activeSet
  idx = .Call("R_retrieve_beta", fit$coefficients, groupSizes, totalGroups, integer(totalGroups), integer(length(fit$coefficients)))
  beta = c(fit$mu, fit$coefficients[idx$betaIdx != 0])
  range = c(0, cumsum(numGroups))
  activeSet = lapply(1:5, function(i){
    if (numGroups[i] > 0){
      index = which(idx$idx[(range[i]+1):range[i+1]] != 0)
      if (length(index) > 0) matrix(activeSet[[i]][index, ], nrow=length(index))
      else NULL
    }
    else NULL
  })
  names(activeSet) = c("cat", "cont", "catcat", "contcont", "catcont")

  #output
  list(betahat=beta, activeSet=activeSet, res=res, objValue=objValue,
       converged=as.logical(fit$converged), iterations=as.integer(fit$iterations))
}




