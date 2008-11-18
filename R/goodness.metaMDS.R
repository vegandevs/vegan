"goodness.metaMDS" <-
    function(object, dis, ...)
{
    require(MASS) || stop("Needs MASS package")
    if (missing(dis))
        dis <- metaMDSredist(object)
    if(attr(dis, "Size") != nrow(object$points))
        stop("Dimensions do not match in ordination and dissimilarities")
    d <- order(dis)
    shep <- Shepard(dis, object$points)
    res <- (shep$y - shep$yf)^2/sum(shep$y^2)
    stress <- sqrt(sum(res))*100
    if ( abs(stress - object$stress) > 0.001)
        stop("Dissimilarities and ordination do not match")
    res <- res[order(d)]
    attr(res, "Size") <- attr(dis, "Size")
    attr(res, "Labels") <- attr(dis, "Labels")
    class(res) <- "dist"
    sqrt(colSums(as.matrix(res))/2*10000)
}

