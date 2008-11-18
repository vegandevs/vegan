"print.procrustes" <-
  function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Procrustes sum of squares:\n")
  cat(formatC(x$ss, digits = digits), "\n\n")
  invisible(x)
}
