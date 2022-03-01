`plot.spantree` <-
    function (x, ord, cex = 0.7, type = "p", labels, dlim, FUN = sammon, 
              ...) 
{
    FUNname <- deparse(substitute(FUN))
    FUN <- match.fun(FUN)
    n <- x$n
    if (missing(ord)) {
        d <- cophenetic(x)
        if (any(d<=0))
            d[d<=0] <- min(d>0)/10
        if (!missing(dlim)) 
            d[d > dlim ] <- dlim
        if (n > 2) {
            ## sammon needs extra care, for other cases we just try FUN(d)
            if (FUNname == "sammon") {
                y <- cmdscale(d)
                dup <- duplicated(y)
                if (any(dup))
                    y[dup, ] <- y[dup,] + runif(2*sum(dup), -0.01, 0.01)
                ord <- FUN(d, y = y)
            } else
                ord <- FUN(d)
        } else
            ord <- cbind(seq_len(n), rep(0,n))
    }
    ord <- scores(ord, display = "sites", ...)
    ordiArgAbsorber(ord, asp = 1, type = "n", FUN = "plot", ...)
    lines(x, ord, ...)
    if (type == "p" || type == "b") 
        ordiArgAbsorber(ord, cex = cex, FUN = "points", ...)
    else if (type == "t") {
        if (missing(labels)) 
            labels <- x$labels
        x <- scores(ord, display = "sites", ...)
        ordiArgAbsorber(x, labels = labels, cex = cex, FUN = "ordilabel", ...)
    }
    ord <- list(sites = ord)
    class(ord) <- "ordiplot"
    invisible(ord)
}
