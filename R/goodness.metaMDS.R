`goodness.metaMDS` <-
    function(object, dis, ...)
{
    if (inherits(object, "monoMDS"))
        return(NextMethod("goodness", object, ...))
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

`goodness.monoMDS` <-
    function(object, ...)
{
    ## Return vector 'x' for which sum(x^2) == stress
    stresscomp <- function(y, yf, form)
    {
        num <- (y-yf)^2
        if (form == 1)
            den <- sum(y^2)
        else
            den <- sum((y-mean(y))^2)
        num/den
    }
    ## Global, local
    if (object$model %in% c("global", "linear")) {
        x <- stresscomp(object$dist, object$dhat, object$isform)
        mat <- matrix(0, object$nobj, object$nobj)
        for (i in 1:object$ndis)
            mat[object$iidx[i], object$jidx[i]] <- x[i]
        res <- sqrt(colSums(mat + t(mat))/2)
    }
    ## Local: returns pointwise components of stress
    else if (object$model == "local") {
        res <- object$grstress/sqrt(object$ngrp)
    } else if (object$model == "hybrid" && object$ngrp == 2) {
        mat <- matrix(0, object$nobj, object$nobj)
        gr <- seq_len(object$ndis) < object$istart[2]
        x <- stresscomp(object$dist[gr], object$dhat[gr], object$isform)
        x <- c(x, stresscomp(object$dist[!gr], object$dhat[!gr], object$isform))
        i <- object$iidx
        j <- object$jidx
        for (k in 1:object$ndis) {
            mat[i[k], j[k]] <- mat[i[k], j[k]] + x[k]
        }
        res <- sqrt(colSums(mat + t(mat))/4)
    } else {
        stop("unknown 'monoMDS' model")
    }
    res
}
