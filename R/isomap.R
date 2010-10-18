`isomap` <-
function(dist, ndim=10, ...)
{
    dist <- isomapdist(dist, ...)
    out <- cmdscale(dist, k=ndim, eig=TRUE)
    ## some versions of cmdscale may return NaN points corresponding
    ## to negative eigenvalues
    if (any(out$eig < 0)) {
        out$eig <- out$eig[out$eig > 0]
        out$points <- out$points[, seq_along(out$eig), drop = FALSE]
        warning(gettextf("isomap returns only %d axes with positive eigenvalues",
                         length(out$eig)))
    }
    npoints <- nrow(out$points)
    net <- matrix(FALSE, nrow=npoints, ncol=npoints)
    net[lower.tri(net)][attr(dist, "net")] <- TRUE 
    net <- which(net, arr.ind=TRUE)
    out$method <- attr(dist, "method")
    out$criterion <- attr(dist, "criterion")
    out$critval <- attr(dist, "critval")
    out$take <- attr(dist, "take")
    out$net <- net
    out$npoints <- npoints
    out$call <- match.call()
    class(out) <- "isomap"
    out
}

