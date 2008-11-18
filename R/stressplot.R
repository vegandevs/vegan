"stressplot" <-
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
        if (length(dis) > 5000) pch = "." else pch = 1
    plot(shep, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    lines(shep$x, shep$yf, type = "S", col = l.col, lwd = lwd, ...)
    lab <- paste("Non-metric fit, R2 =", format(rstress, digits=3),
               "\nLinear fit, R2 =", format(ralscal, digits=3))
    text(min(shep$x), 0.95*max(shep$y), lab, pos=4)
    invisible(shep)
}
