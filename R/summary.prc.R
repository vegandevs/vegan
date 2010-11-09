`summary.prc` <-
    function (object, axis = 1, scaling = 3, digits = 4, ...)
{
    sc = scores(object, scaling = scaling, display = c("sp", "lc"),
                           choices=axis, ...)
    ## coef for scaled sites (coef(object) gives for orthonormal)
    b <- qr.coef(object$CCA$QR, sc$constraints)
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
