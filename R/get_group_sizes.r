get_group_sizes = function(activeSet, numLevels){

  if (is.null(activeSet$cat)) numCat = NULL
  else numCat = numLevels[activeSet$cat]

  if (is.null(activeSet$cont)) numCont = NULL
  else numCont = rep(1, length(activeSet$cont))
  
  if (is.null(activeSet$catcat)) numCatCat = NULL
  else numCatCat = apply(activeSet$catcat, 1, function(x) prod(numLevels[x]))
  
  if (is.null(activeSet$contcont)) numContCont = NULL
  else numContCont = rep(3, nrow(activeSet$contcont))
  
  if (is.null(activeSet$catcont)) numCatCont = NULL
  else numCatCont = 2*numLevels[activeSet$catcont[, 1]]
  
  sizes = c(numCat, numCont, numCatCat, numContCont, numCatCont)
  if (any(!is.finite(sizes)) || any(sizes > .Machine$integer.max) ||
      sum(as.double(sizes)) + 1 > .Machine$integer.max)
    stop("active coefficient vector exceeds native integer limits")
  as.integer(sizes)
}
