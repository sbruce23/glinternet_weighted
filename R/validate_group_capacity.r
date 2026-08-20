validate_group_capacity = function(numLevels, interactionCandidates=NULL, interactionPairs=NULL) {
  intMax = .Machine$integer.max
  categorical = which(numLevels > 1)
  if (sum(as.double(numLevels[categorical])) + 1 > intMax)
    stop("declared categorical main-effect coefficient length exceeds native integer limits")

  checkCatCat = function(left, right) {
    if (numLevels[left] > intMax / numLevels[right])
      stop("a categorical interaction group exceeds native integer limits")
  }
  checkCatCont = function(cat) {
    if (numLevels[cat] > intMax / 2)
      stop("a categorical-continuous interaction group exceeds native integer limits")
  }

  if (!is.null(interactionPairs)) {
    for (i in seq_len(nrow(interactionPairs))) {
      left = interactionPairs[i,1]
      right = interactionPairs[i,2]
      if (numLevels[left] > 1 && numLevels[right] > 1) checkCatCat(left,right)
      else if (numLevels[left] > 1) checkCatCont(left)
      else if (numLevels[right] > 1) checkCatCont(right)
    }
    return(invisible(NULL))
  }

  eligible = if (is.null(interactionCandidates)) seq_along(numLevels) else interactionCandidates
  eligibleCat = eligible[numLevels[eligible] > 1]
  eligibleCont = eligible[numLevels[eligible] == 1]
  for (cat in eligibleCat) {
    partners = setdiff(categorical, cat)
    if (length(partners)) checkCatCat(cat, partners[which.max(numLevels[partners])])
  }
  continuous = which(numLevels == 1)
  if (length(eligibleCat) && length(continuous))
    checkCatCont(eligibleCat[which.max(numLevels[eligibleCat])])
  if (length(eligibleCont) && length(categorical))
    checkCatCont(categorical[which.max(numLevels[categorical])])
  invisible(NULL)
}
