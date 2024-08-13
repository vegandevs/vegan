`ca` <- function(X, ...) {
  if (inherits(X, "formula")) {
    stop("Argument 'X' was supplied a formula, which is not supported by 'ca()'",
    call. = FALSE)
  }
  ord <- cca(X = X, ...)
  # change the call to be from ca()
  ord$call <- match.call()
  class(ord) <- append(class(ord), "vegan_ca", after = 0)
  ord
}
