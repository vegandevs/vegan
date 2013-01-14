`print.mantel` <-
  function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\n")
  if (inherits(x, "mantel.partial"))
    cat("Partial ")
  cat("Mantel statistic based on", x$method, "\n")
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Mantel statistic r: ")
  cat(formatC(x$statistic, digits = digits), "\n")
  nperm <- x$permutations
  if (nperm) {
    cat("      Significance:", format.pval(x$signif), 
        "\n\n")
    out <- quantile(x$perm, c(0.9, 0.95, 0.975, 0.99))
    cat("Upper quantiles of permutations (null model):\n")
    print(out, digits = 3)
    cat("\nBased on", nperm, "permutations")
    if (!is.null(x$strata)) 
      cat(", stratified within", x$strata)
  }
  cat("\n\n")
  invisible(x)
}

