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
  invisible(x)
}

