`spenvcor` <-
    function (object)
{
    if (is.null(object$CCA))
        stop("Needs results from constrained ordination")
    u <- object$CCA$u
    wa <- object$CCA$wa
    if (!inherits(object, "rda")) { # is CCA
        r <- sqrt(object$rowsum)
        u <- r * u
        wa <- r * wa
    }
    ## because colSums(u*u) = 1, we can simplify diag(cor(u, wa)) --
    ## and we must for weighted CA
    colSums(u * wa)/sqrt(colSums(wa^2))
}

