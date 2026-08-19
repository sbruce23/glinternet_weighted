rescale_betahat = function(activeSet, betahat, X, Z, weights, levels, n){

  if (is.null(activeSet)) return(betahat)
  
  nVars = sapply(activeSet, function(x) if (is.null(x)) 0 else nrow(x))
  betaLen = length(betahat)
  indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)

  result = .Call("R_rescale_beta", X, Z, weights, n, betahat, betaLen, nVars, levels, indices$cat, indices$cont, indices$catcat, indices$contcont, indices$catcont, double(betaLen))
  if (nVars["contcont"] > 0) {
    sizes = get_group_sizes(activeSet, levels)
    starts = 2 + c(0,cumsum(sizes)[-length(sizes)])
    firstContCont = nVars["cat"] + nVars["cont"] + nVars["catcat"] + 1
    groups = seq.int(firstContCont, length.out=nVars["contcont"])
    interactionPositions = starts[groups] + 2
    unrepresentable = betahat[interactionPositions] != 0 &
      (!is.finite(result[interactionPositions]) | result[interactionPositions] == 0)
    if (any(unrepresentable))
      stop("a fitted continuous interaction coefficient is not representable on the original predictor scale")
  }
  result
}

  
