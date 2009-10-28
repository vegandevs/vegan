"pcnm" <-
    function(matdist, threshold, support = c("vegan", "ade4"), w)
{
    EPS <- sqrt(.Machine$double.eps)
    wa.old <- options(warn = -1)
    on.exit(options(wa.old))
    matdist <- as.dist(matdist)
    if (missing(threshold)) {
        support <- match.arg(support)
        if (!missing(w) && support == "ade4")
            stop("weights are not supported with 'ade4'")
        threshold <- 
            switch(support,
                   vegan =  max(spantree(matdist)$dist),
                   ade4 = max(neig2mat(mstree(matdist)) * as.matrix(matdist))
                   )
    }
    matdist[matdist > threshold] <- 4*threshold
    ## vegan:::wcmdscale used to be able to use weights which also
    ## means that 'k' need not be given, but all vecctorw with >0
    ## eigenvalues will be found
    mypcnm <- wcmdscale(matdist, eig = TRUE, w=w)
    res <- list(vectors = mypcnm$points, values = mypcnm$eig,
                weights = mypcnm$weig)
    k <- ncol(mypcnm$points)
    res$vectors <- sweep(res$vectors, 2, sqrt(res$values[1:k]), "/")
    colnames(res$vectors) <- paste("PCNM", 1:k, sep="")
    res$threshold <- threshold
    class(res) <- "pcnm"
    res
}
