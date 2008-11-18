`isomap` <-
function(dist, ndim=10, ...)
{
    dist <- isomapdist(dist, ...)
    out <- cmdscale(dist, k=ndim, eig=TRUE)
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

