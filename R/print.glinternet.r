print.glinternet = function(x, ...){

  groupNames = c("cat", "cont", "catcat", "contcont", "catcont")
  sizes = t(vapply(x$activeSet, function(active) {
    if (is.null(active)) return(setNames(integer(5), groupNames))
    vapply(groupNames, function(name) {
      group = active[[name]]
      if (is.null(group)) 0L else nrow(group)
    }, integer(1))
  }, setNames(integer(5), groupNames)))
  output = data.frame(lambda=signif(x$lambda, 3), objValue=signif(x$objValue, 3), sizes)

  #print the call
  cat("Call: ")
  print(x$call)

  #print number of interactions found
  print(output)
}
