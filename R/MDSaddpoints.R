### Add New Points to NMDS ordination

`MDSaddpoints` <-
    function (nmds, dis, neighbours=5, maxit=200)
{
## bring list component to local matrix

    points <- nmds$points
    oldn <- nrow(points)
    ndim <- ncol(points)

    ## use dist2xy to handle input 'dis'
    if (inherits(dis, "dist")) {
        totn <- attr(dis, "Size")
        if (totn <= oldn)
            stop("input dissimilarities should have more observations than 'nmds'")
        dis <- dist2xy(dis, pick = seq(oldn + 1, totn))
    }
    newn <- nrow(dis)
    totn <- oldn + newn

    ## set up initial coordinates as weighted average of nearest
    ## neighbours to old points using 1-diss as weights

    tmp <- matrix(0, newn, ndim)
    for (i in 1:newn) {
        pnt <- order(dis[i,seq_len(oldn)])[seq_len(neighbours)]
        maxdist <- attr(dis, "maxdist")
        if (is.null(maxdist))
            maxdist <- max(1, dis)
        weight <- maxdist - dis[i,pnt]
        for (j in 1:ncol(points)) {
            tmp[i,j] <- weighted.mean(points[pnt,j], w=weight)
        }
    }

    xinit <- rbind(points,tmp)

    ## set up indices

    iidx <- as.vector(row(dis)) + oldn
    jidx <- as.vector(col(dis))

    ## combine with old data

    diss <- c(nmds$diss, dis)
    ndis <- length(diss)
    iidx <- c(nmds$iidx, iidx)
    jidx <- c(nmds$jidx, jidx)

    ## set up ordination parameters.
    nfix <- oldn
    ngrp <-
    istart <- 1
    dist <- rep(0,ndis)
    dhat <- rep(0,ndis)
    x <- matrix(0,nrow=totn,ncol=ndim)
    stress <- 1.0
    strs <- ngrp
    iters <- 1
    icause <- 1
    maxits <- maxit
    iwork <- rep(0,ndis)
    grad <- matrix(0,nrow=totn,ncol=ndim)
    grlast <- matrix(0,nrow=totn,ncol=ndim)

    out <- .Fortran('monoMDS',
           nobj=as.integer(totn),
           nfix=as.integer(oldn),
           ndim=as.integer(ndim),
           ndis=as.integer(ndis),
           ngrp=as.integer(1),
           diss=as.double(diss),
           iidx=as.integer(iidx),
           jidx=as.integer(jidx),
           xinit=as.double(xinit),
           istart=as.integer(1),
           isform=as.integer(nmds$isform),
           ities=as.integer(nmds$ities),
           iregn=as.integer(nmds$iregn),
           iscal=as.integer(0), # should be 0 with nfix
           maxits=as.integer(maxits),
           sratmx=as.double(nmds$sratmx),
           strmin=as.double(nmds$strmin),
           sfgrmn=as.double(nmds$sfgrmn),
           dist=as.double(dist),
           dhat=as.double(dhat),
           points=as.double(x),
           stress=as.double(stress),
           strs=as.double(1),
           iters=as.integer(iters),
           cause=as.integer(icause))
    dim(out$points) <- c(totn, ndim)
    newpoints <- out$points[(oldn+1):totn,, drop=FALSE]
    dimnames(newpoints) <- list(rownames(dis), colnames(nmds$points))
    adds <- list(points = newpoints, seed = tmp,
                 deltastress = out$stress - nmds$stress,
                 iters = out$iters, cause = out$cause)
    class(adds) <- 'nmds'
    adds
}
