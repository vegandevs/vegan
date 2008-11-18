"summary.prc" <-
    function (object, axis = 1, scaling = 2, digits = 4, ...) 
{
    species <- drop(scores(object, scaling = scaling, display="sp", choices=axis))
    b <- coef(object)[, axis]
    prnk <- object$pCCA$rank
    lentreat <- length(object$terminfo$xlev[[2]])
    lenb <- length(b)
    b <- b[-(1:(2 * prnk))]
    bx <- b[1:(lentreat - 1)]
    by <- b[lentreat:length(b)]
    b <- cbind(bx, matrix(by, nrow = lentreat - 1, byrow = TRUE))
    rownames(b) <- (object$terminfo$xlev[[2]])[-1]
    colnames(b) <- object$terminfo$xlev[[1]]
    out <- list(sp = species, coefficients = b, names = names(object$terminfo$xlev), 
                corner = (object$terminfo$xlev[[2]])[1], call = object$call, 
                digits = digits)
    class(out) <- "summary.prc"
    out
}
