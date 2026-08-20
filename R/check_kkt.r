check_kkt = function(X, Z, res, weights, n, pCat, pCont, numLevels, candidates, activeSet, lambda, numCores=1){

  norms = vector("list", 5)
  names(norms) = c("cat", "cont", "catcat", "contcont", "catcont")

  #compute the norms
  if (pCat > 0){
    norms$cat = compute_norms_cat(X, weights * res, n, pCat, numLevels, numCores)
    if (!is.null(candidates$variables$catcat)) norms$catcat = compute_norms_cat_cat(X, weights * res, n, numLevels, candidates$variables$catcat, numCores)
  }
  if (pCont > 0){
    norms$cont = compute_norms_cont(Z, weights * res, n)
    if (!is.null(candidates$variables$contcont)) norms$contcont = compute_norms_cont_cont(Z, norms$cont, weights * res, weights, n, candidates$variables$contcont, numCores)
  }
  if (!is.null(candidates$variables$catcont)) norms$catcont = compute_norms_cat_cont(X, Z, norms$cat, weights * res, n, numLevels, candidates$variables$catcont, numCores)

  #check for nonzero variables
  violators = lapply(names(norms), function(x) {
    if (!is.null(norms[[x]])) {
      indices = which(norms[[x]] > lambda)
      if (length(indices) > 0) {
        return (matrix(candidates$variables[[x]][indices, ], nrow=length(indices)))
      }
    }
    return (NULL)
  })
  names(violators) = names(norms)

  #check if any variables should be added to active set
  flag = 1
  for (nm in names(norms)) {
    if (is.null(violators[[nm]])) {
      next
    }
    if (!is.null(activeSet[[nm]])) {
      actives = apply(activeSet[[nm]], 1, function(x) paste(x, collapse=":"))
      extras = which(!(apply(violators[[nm]], 1, function(x) paste(x, collapse=":")) %in% actives))
    } else {
      extras = 1:nrow(violators[[nm]])
    }
    if (length(extras) > 0) {
      # Preserve a matrix when multiple rows are selected from a one-column
      # main-effect group. Without drop=FALSE, R turns those rows into a vector;
      # rbind then recycles/truncates it and can silently omit KKT violators.
      activeSet[[nm]] = rbind(
        activeSet[[nm]],
        violators[[nm]][extras, , drop=FALSE]
      )
      flag = 0
    }
  }

  list(norms=norms, activeSet=activeSet, flag=flag)
}
    
    
