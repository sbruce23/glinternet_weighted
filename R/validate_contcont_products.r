validate_contcont_products = function(X, pairs, weights=NULL) {
  if (is.null(pairs) || !nrow(pairs)) return(invisible(NULL))
  rows = if (is.null(weights)) rep(TRUE,nrow(X)) else weights > 0
  limit = .Machine$double.xmax
  for (i in seq_len(nrow(pairs))) {
    left = abs(X[rows,pairs[i,1]])
    right = abs(X[rows,pairs[i,2]])
    nonzero = right > 0
    if (any(nonzero & left > limit/right))
      stop("continuous interaction contains a raw product outside the finite double range")
  }
  invisible(NULL)
}
