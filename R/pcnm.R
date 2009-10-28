"pcnm" <-
    function(matdist, threshold, support = c("vegan", "ade4"))
{
    EPS <- sqrt(.Machine$double.eps)
    wa.old <- options(warn = -1)
    on.exit(options(wa.old))
    matdist <- as.dist(matdist)
    if (missing(threshold)) {
        support <- match.arg(support)
        threshold <- 
            switch(support,
                   vegan =  max(spantree(matdist)$dist),
                   ade4 = max(neig2mat(mstree(matdist)) * as.matrix(matdist))
                   )
    }
    matdist[matdist > threshold] <- 4*threshold
    k <- attr(matdist, "Size") - 1
    mypcnm <- cmdscale(matdist, k = k, eig = TRUE)
    eq0 <- abs(mypcnm$eig/max((mypcnm$eig))) <= EPS
    inf0 <- mypcnm$eig < 0
    res <- list()
    res$values <- mypcnm$eig[!(eq0|inf0)]
    res$vectors <- mypcnm$points[,!(eq0|inf0), drop = FALSE]
    res$vectors <- sweep(res$vectors, 2, sqrt(res$values), "/")
    colnames(res$vectors) <- paste("PCNM", 1:ncol(res$vectors), sep="")
    res$threshold <- threshold
    class(res) <- "pcnm"
    res
}
