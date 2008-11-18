"plot.spantree" <-
    function (x, ord, cex = 0.7, type = "p", labels, dlim, FUN = sammon, 
              ...) 
{
    FUNname <- deparse(substitute(FUN))
    if (length(FUNname) && FUNname %in% c("sammon", "isoMDS")) 
        require(MASS) || stop(FUNname, "requires package MASS")
    FUN <- match.fun(FUN)
    n <- length(x$kid) + 1
    if (missing(ord)) {
        d <- cophenetic(x)
        if (any(d<=0))
            d[d<=0] <- min(d>0)/10
        if (!missing(dlim)) 
            d[d > dlim ] <- dlim
        y <- cmdscale(d)
        dup <- duplicated(y)
        if (any(dup))
            y[dup, ] <- y[dup,] + runif(2*sum(dup), -0.01, 0.01) 
        ord <- FUN(d, y)
    }
    ord <- scores(ord, display = "sites")
    plot(ord, asp = 1, type = "n", ...)
    lines(x, ord)
    if (type == "p" || type == "b") 
        points(ord, cex = cex, ...)
    else if (type == "t") {
        if (missing(labels)) 
            labels <- x$labels
        ordilabel(ord, display = "sites", labels = labels, cex = cex, ...)
    }
    ord <- list(sites = ord)
    class(ord) <- "ordiplot"
    invisible(ord)
}
