`pco` <- function(X, ...) {
  if (inherits(X, "formula")) {
    stop("Argument 'X' was supplied a formula, which is not supported by 'pco()'",
    call. = FALSE)
  }
  ord <- dbrda(X ~ 1, ...)
  # change the call to be from pco()
  ord$call <- match.call()
  class(ord) <- append(class(ord), "vegan_pco", after = 0)
  ord
}
