addpoints.nmds <- function (nmds,dis,neighbors=5,maxits=200)
{
# bring list component to local matrix

    points <- nmds$points

# get sizes

    oldn <- nrow(points)
    ndim <- ncol(points)
    totn <- attr(dis,'Size')
    newn <- totn - oldn

# test correspondence

    if (!identical(dimnames(points)[[1]],attr(dis,'Labels')[1:oldn]))
         stop('ordination and dissimilarity matrix do not match')

# decompose disimilarity object to isolate new values

    diss <- as.matrix(dis)[1:oldn,(oldn+1):totn]
    
# set up initial conditions

    ndis <- oldn * newn
    tmp <- matrix(rep(0,newn*ndim),ncol=ndim)
    for (i in 1:newn) {
        pnt <- seq(1,oldn)[order(diss[,i])][1:neighbors]
        weight <- 1-diss[pnt,i]
        for (j in 1:ncol(points)) {
            tmp[i,j] <- weighted.mean(points[pnt,j],w=weight)
        }
    }
    
    xinit <- rbind(points,tmp)
    dimnames(xinit)[[1]] <- attr(dis,'Labels')

# set up indices

    iidx <- rep((1:oldn),newn)
    jidx <- NULL
    for (i in (oldn+1):totn) jidx <- c(jidx,rep(i,oldn))

#set up ordination

    nfix <- oldn
    ngrp <- 
    istart <- 1
    isform <- 2
    ities <- 1
    iregn <- 1
    iscal <- 0
    sratmx <- 0.99999
    strmin <- 1e-07
    sfgrmn <- 1e-05
    dist <- rep(0,ndis)
    dhat <- rep(0,ndis)
    x <- matrix(0,nrow=totn,ncol=ndim)
    stress <- 1
    strs <- ngrp
    iters <- 1
    icause <- 1 
    maxits <- maxits
    iwork <- rep(0,ndis)
    grad <- matrix(0,nrow=totn,ncol=ndim)
    grlast <- matrix(0,nrow=totn,ncol=ndim)

    out <- .Fortran('monoMDS',
           nobj=as.integer(totn),
           nfix=as.integer(nfix),
           ndim=as.integer(ndim),
           ndis=as.integer(ndis),
           ngrp=as.integer(ngrp),
           diss=as.double(diss),
           iidx=as.integer(iidx),
           jidx=as.integer(jidx),
           xinit=as.double(xinit),
           istart=as.integer(istart),
           isform=as.integer(isform),
           ities=as.integer(ities),
           iregn=as.integer(iregn),
           iscal=as.integer(iscal),
           maxits=as.integer(maxits),
           sratmx=as.double(sratmx),
           strmin=as.double(strmin),
           sfgrmn=as.double(sfgrmn),
           dist=as.double(dist),
           dhat=as.double(dhat),
           points=as.double(x),
           stress=as.double(stress),
           strs=as.double(strs),
           iters=as.integer(iters),
           cause=as.integer(icause))
    res <- list(points=matrix(out$points,ncol=ndim),
                newpoints=matrix(out$points,ncol=ndim)[(oldn+1):totn,],
                stress=out$stress,iters=out$iters,cause=out$cause,
                seeds=tmp)
    dimnames(res$points)[[1]] <- attr(dis,'Labels')       
    dimnames(res$newpoints)[[1]] <- attr(dis,'Labels')[(oldn+1):totn]
    class(res) <- 'nmds'
    res 
}
