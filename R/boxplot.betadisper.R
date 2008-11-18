`boxplot.betadisper` <- function(x, ylab = "Distance to centroid", ...) {
  tmp <- boxplot(x$distances ~ x$group, ylab = ylab, ...)
  invisible(tmp)
}
