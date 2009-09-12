"print.envfit" <-
  function(x, ...)
{
  if (!is.null(x$vectors)) {
    cat("\n***VECTORS\n\n")
    print(x$vectors)
  }
  if (!is.null(x$factors)) {
    cat("\n***FACTORS:\n\n")
    print(x$factors)
  }
  if (!is.null(x$na.action))
      cat("\n", naprint(x$na.action), "\n", sep="")
  invisible(x)
}

