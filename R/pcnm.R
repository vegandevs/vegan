"pcnm" <-
function(matdist,thresh=give.thresh(as.dist(matdist)))
{
    matdist <- as.matrix(matdist)
   
   
            mattrunc <- ifelse(matdist >thresh, 4*thresh,matdist)
    wa.old <- options()$warn
    options(warn = -1)
    mypcnm <- cmdscale(mattrunc,k=min(dim(matdist))-1,eig=TRUE)
    eq0 <- apply(as.matrix(mypcnm$eig/max((mypcnm$eig))),1,function(x) identical(all.equal(x, 0), TRUE))
    inf0 <- ifelse(mypcnm$eig<0,TRUE,FALSE)
    res <- list()
    res$values <- mypcnm$eig[!(eq0|inf0)]
    res$vectors <- mypcnm$points[,!(eq0|inf0)]
    res$vectors <- sweep(res$vectors,2,sqrt(res$values),"/")
    options(warn = wa.old)
    return(res)
}
