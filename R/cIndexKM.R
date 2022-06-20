`cIndexKM` <-
    function (y, x, index = "all")
{
    kmeans_res <- y
#############################################
    gss <- function(x, clsize, withins)
    {
        allmean <- colMeans(x)
        dmean <- sweep(x, 2, allmean, "-")
        allmeandist <- sum(dmean^2)
        wgss <- sum(withins)
        bgss <- allmeandist - wgss
        list(wgss = wgss, bgss = bgss)
    }
#############################################

### Function modified by SD and PL from the original "cIndexKM" in "cclust"
### to accommodate a single response variable as well as singleton groups
### and remove unwanted index.

###  The index
################################################
    calinski <- function(zgss, clsize)
    {
        n <- sum(clsize)
        k <- length(clsize)
        ## undefined 0/0 for one class (or fewer in error cases)
        if (k <= 1)
            NA
        else
            zgss$bgss/(k - 1)/(zgss$wgss/(n - k))
    }
################################################
    ssi <- function(centers, clsize)
    {
        ncl <- dim(centers)[1]
        nvar <- dim(centers)[2]
        cmax <- apply(centers, 2, max)
        cmin <- apply(centers, 2, min)
        cord <- apply(centers, 2, order)
        cmaxi <- cord[ncl, ]
        cmini <- cord[1, ]
        meanmean <- mean(centers)
        absmdif <- abs(apply(centers, 2, mean) - meanmean)
        span <- cmax - cmin
        csizemax <- clsize[cmaxi]
        csizemin <- clsize[cmini]
        hiest <- nvar
        hiestw <- hiest * max(max(csizemax), max(csizemin)) *
            exp(-min(absmdif))
        sist <- sum(span)/hiest
        sistw <- (span * exp(-absmdif)) %*% sqrt(csizemax * csizemin)/hiestw
        list(ssi = sist, ssiw = sistw)
    }
################################################

    zgss <- gss(x, kmeans_res$size, kmeans_res$withinss)

    index <- pmatch(index, c("calinski", "ssi", "all"))
    if (is.na(index))
        stop("invalid clustering index")
    if (index == -1)
        stop("ambiguous index")
    vecallindex <- numeric(3)
    if (any(index == 1) || (index == 3))
        vecallindex[1] <- calinski(zgss, kmeans_res$size)
    if (any(index == 2) || (index == 3))
        vecallindex[2] <- ssi(kmeans_res$centers, kmeans_res$size)$ssiw
    names(vecallindex) <- c("calinski", "ssi")
    if (index < 3)
        vecallindex <- vecallindex[index]
    vecallindex
}
