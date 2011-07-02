`stressplot`<-
    function(object, ...)
{
    UseMethod("stressplot")
}

`stressplot.monoMDS` <-
    function(object, pch, p.col = "blue", l.col = "red", lwd, ...)
{
    if (missing(lwd))
        if (object$ngrp > 2)
            lwd <- 1
        else
            lwd <- 2
    ## extract items to plot
    x <- object$diss
    y <- object$dist
    yf <- object$dhat
    ## all models plot dist against diss, but there can be duplicated
    ## items in some models: remove duplicates in hybrid (iregn==3)
    ## and local (ngrp > 1) models:
    if (object$iregn == 3)
        pts <- seq_along(x) < object$istart[2]
    else if (object$ngrp > 2)
        pts <- object$iidx > object$jidx
    else
        pts <- !logical(length(x))
    ## Plotting character
    if (missing(pch))
        if (sum(pts) > 5000) pch <- "." else pch <- 1
    ## plot points
    plot(x[pts], y[pts], pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    ## collect values for 'linear fit'
    ralscal <- 0
    ## Fit lines: linear (iregn=2) and hybrid (iregn=3) have a smooth line
    if (object$iregn > 1) {
        if (object$iregn == 3) {
            k <- seq(object$istart[2], object$ndis)
            yl <- range(yf[k])
            xl <- range(x[k])
            ralscal <- cor(y[k], yf[k])^2
        } else {
            yl <- range(yf)
            xl <- range(x)
            ralscal <- cor(y, yf)^2
        }
        lines(xl, yl, col = l.col, lwd = lwd, ...)
    }
    ## Monotone line except in linear, and local has several...
    if (object$iregn != 2) {
        ist <- c(object$istart, object$ndis + 1)
        if (object$iregn == 3)
            object$ngrp <- 1
        for(j in 1:object$ngrp) {
            k <- seq(ist[j], ist[j+1]-1)
            ralscal <- ralscal + cor(y[k], yf[k])^2
            lines(x[k], yf[k], type = "S", col = l.col, lwd = lwd, ...)
        }
    }
    ## Stress as R2
    rstress <- 1 - object$stress^2
    ralscal <- if(object$iregn == 3) ralscal/2 else ralscal/object$ngrp
    lab <- paste("Non-metric fit, R2 =", format(rstress, digits=3),
                 "\nLinear fit, R2 =", format(ralscal, digits=3))
    text(min(x), 0.95*max(y), lab, pos=4)
    invisible(list("x" = x, "y" = y, "yf" = yf))
}
    
`stressplot.default` <-
    function(object, dis, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    require(MASS) || stop("Needs MASS package")
    if (missing(dis))
        dis <- metaMDSredist(object)
    if (attr(dis, "Size") != nrow(object$points))
        stop("Dimensions do not match in ordination and dissimilarities")
    shep <- Shepard(dis, object$points)
    stress <- sum((shep$y - shep$yf)^2)/sum(shep$y^2)
    rstress <- 1 - stress
    ralscal <- cor(shep$y, shep$yf)^2
    stress <- sqrt(stress)*100
    if ( abs(stress - object$stress) > 0.001)
        stop("Dissimilarities and ordination do not match")
    if (missing(pch))
        if (length(dis) > 5000) pch <- "." else pch <- 1
    plot(shep, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    lines(shep$x, shep$yf, type = "S", col = l.col, lwd = lwd, ...)
    lab <- paste("Non-metric fit, R2 =", format(rstress, digits=3),
               "\nLinear fit, R2 =", format(ralscal, digits=3))
    text(min(shep$x), 0.95*max(shep$y), lab, pos=4)
    invisible(shep)
}
