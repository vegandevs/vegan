"cIndexKM" <- function (y, x, index = "all") 
{
    kmeans_res <- y
#########################################
    withinss <- function(kmeans_res, x) 
    {
        retval <- rep(0, nrow(kmeans_res$centers))
        x <- (x - kmeans_res$centers[kmeans_res$cluster, ])^2
        for (k in 1:nrow(kmeans_res$centers)) 
        {
            retval[k] <- sum(x[kmeans_res$cluster == k, ])
        }
        retval
    }
##########################################
    varwithinss <- function(x, centers, cluster) 
    {
        nrow <- dim(centers)[1]
        nvar <- dim(x)[2]
        varwithins <- matrix(0, nrow, nvar)
        x <- (x - centers[cluster, ])^2
        for (l in 1:nvar) 
        {
            for (k in 1:nrow) 
            {
                varwithins[k, l] <- sum(x[cluster == k, l])
            }
        }
        return(varwithins)
    }
##########################################
    maxmindist <- function(clsize, distscen) 
    {
        ncl <- length(clsize)
        npairs <- 0
        for (i in 1:ncl) npairs <- npairs + clsize[i] * (clsize[i] - 1)/2
        mindw <- 0
        nfound <- distscen[1]
        i <- 1
        while (nfound < npairs) 
        {
            if ((nfound + distscen[i + 1]) < npairs) 
            {
                mindw <- mindw + i * distscen[i + 1]
                nfound <- nfound + distscen[i + 1]
            }
            else 
            {
                mindw <- mindw + i * (npairs - nfound)
                nfound <- nfound + distscen[i + 1]
            }
            i <- i + 1
        }
        maxdw <- 0
        nfound <- 0
        i <- length(distscen) - 1
        while (nfound < npairs) 
        {
            if ((nfound + distscen[i + 1]) < npairs) 
            {
                maxdw <- maxdw + i * distscen[i + 1]
                nfound <- nfound + distscen[i + 1]
            }
            else 
            {
                maxdw <- maxdw + i * (npairs - nfound)
                nfound <- nfound + distscen[i + 1]
            }
            i <- i - 1
        }
        minmaxd <- list(mindw = mindw, maxdw = maxdw)
        return(minmaxd)
    }
#############################################
    gss <- function(x, clsize, withins) 
    {
        n <- sum(clsize)
        k <- length(clsize)
        allmean <- colMeans(x)
        dmean <- sweep(x, 2, allmean, "-")
        allmeandist <- sum(dmean^2)
        wgss <- sum(withins)
        bgss <- allmeandist - wgss
        zgss <- list(wgss = wgss, bgss = bgss)
        return(zgss)
    }
#############################################
    vargss <- function(x, clsize, varwithins) 
    {
        nvar <- dim(x)[2]
        n <- sum(clsize)
        k <- length(clsize)
        varallmean <- rep(0, nvar)
        varallmeandist <- rep(0, nvar)
        varwgss <- rep(0, nvar)
        for (l in 1:nvar) varallmean[l] <- mean(x[, l])
        vardmean <- sweep(x, 2, varallmean, "-")
        for (l in 1:nvar) 
        {
            varallmeandist[l] <- sum((vardmean[, l])^2)
            varwgss[l] <- sum(varwithins[, l])
        }
        varbgss <- varallmeandist - varwgss
        vartss <- varbgss + varwgss
        zvargss <- list(vartss = vartss, varbgss = varbgss)
        return(zvargss)
    }
		
#################################################
    count <- function(x) 
    {
        nr <- nrow(x)
        nc <- ncol(x)
        d <- integer(nc + 1)
        retval <- .C("count", xrows = nr, xcols = nc, x = as.integer(x), 
                     d = d, PACKAGE = "cclust")
        d <- retval$d
        names(d) <- 0:nc
        return(d)
    }
################################################
### Function modified by SD and PL from the original "cIndexKM" in "cclust" 
### to accommodate a single response variable as well as singleton groups
### and remove unwanted index.
		
###  The index
################################################
    calinski <- function(zgss, clsize) 
    {
        n <- sum(clsize)
        k <- length(clsize)
        vrc <- (zgss$bgss/(k - 1))/(zgss$wgss/(n - k))
        return(vrc = vrc)
    }
################################################
    ssi <- function(centers, clsize) 
    {
        ncl <- dim(centers)[1]
        nvar <- dim(centers)[2]
        n <- sum(clsize)
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
        return(list(ssi = sist, ssiw = sistw))
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
