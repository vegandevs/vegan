`pca` <- function(X, scale = FALSE, ...) {
  if (inherits(X, "formula")) {
    stop("Argument 'X' was supplied a formula, which is not supported by 'pca()'",
    call. = FALSE)
  }
  ord <- rda(X = X, scale = scale, ...)
  # change the call to be from pca()
  ord$call <- match.call()
  class(ord) <- append(class(ord), "vegan_pca", after = 0)
  ord
}
