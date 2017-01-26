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
        diag(cor(u, wa)) ## does new centring
    } else { # not CCA, no weights
    ## because sum(u*u) = 1, we can simplify diag(cor(u, wa))
        colSums(u * wa)/sqrt(colSums(wa^2))
    }
}

