`summary.prc` <- function (object, axis = 1, scaling = "symmetric", const,
                           digits = 4, correlation = FALSE, ...) {
    sc = scores(object, scaling = scaling, display = c("sp", "reg"), const,
                choices = axis, correlation = correlation, ...)
    b <- sc$regression
    prnk <- object$pCCA$rank
    lentreat <- length(object$terminfo$xlev[[2]])
    b = matrix(b[-(1:prnk)], nrow = lentreat-1, byrow = TRUE)
    rownames(b) <- (object$terminfo$xlev[[2]])[-1]
    colnames(b) <- object$terminfo$xlev[[1]]
    out <- list(sp = drop(sc$species), coefficients = b,
                names = names(object$terminfo$xlev),
                corner = (object$terminfo$xlev[[2]])[1],
                call = object$call, digits = digits)
    class(out) <- "summary.prc"
    out
}
