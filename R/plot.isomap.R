`plot.isomap` <-
    function (x, net = TRUE, n.col = "gray", type = "points", ...)
{
    type <- match.arg(type, c("points", "text", "none"))
    if (!net) {
        pl <- ordiplot(x, display="sites", type = type, ...)
    } else {
        pl <- ordiplot(x, display = "sites", type = "none", ...)
        z <- scores(pl, display = "sites")
        k <- x$net
        ## recycle colour for points
        n.col <- rep(n.col, length = nrow(z))
        n.col <- col2rgb(n.col)/255
        ## get average of colours of connected points
        n.col <- (n.col[,k[,1]] + n.col[,k[,2]])/2
        n.col <- rgb(t(n.col))
        segments(z[k[,1],1], z[k[,1],2], z[k[,2],1], z[k[,2],2], col=n.col)
        if (type == "points")
            points(pl, "sites", ...)
        else if (type == "text")
            ordilabel(pl, ...)
    }
    invisible(pl)
}

