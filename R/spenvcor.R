`spenvcor` <-
    function (object) 
{
    if (is.null(object$CCA))
        stop("Needs results from constrained ordination")
    u <- object$CCA$u
    wa <- object$CCA$wa
    if (!inherits(object, "rda")) { # is CCA
        r <- sqrt(object$rowsum)
        u <- sweep(u, 1, r, "*")
        wa <- sweep(wa, 1, r, "*")
    }
    diag(cor(u, wa))
}

