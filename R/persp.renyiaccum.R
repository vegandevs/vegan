`persp.renyiaccum` <-
    function(x, theta = 220, col = heat.colors(100), zlim, ...)
{
    dn <- dimnames(x)
    Sites <- seq(along=dn[[1]])
    Scale <- seq(along=dn[[2]])
    Diversity <- x[,,"mean"]
    if (missing(zlim))
        zlim <- range(Diversity, 0)
    if (length(col) > 1) {
        ind <- Diversity
        ind <- (ind[-1,-1] + ind[-1,-ncol(ind)] + ind[-nrow(ind),-1] +
                ind[-nrow(ind),-ncol(ind)])/4
        ind <- round((length(col) - 1) * (ind - min(ind))/diff(range(ind)) + 1)
        col <- col[ind]
    }
    persp(Sites, Scale, Diversity, theta = theta, zlim = zlim, col = col,  ...)
}
