#' Add New Points to NMDS ordination
#'
#' @param nmds Result object from \code{\link{metaMDS}}. The
#'     configuration of points will be kept fixed, but new points will
#'     be added.
#' @param dis Dissimilarity object where the first items are
#'     dissmilarities among fixed points in \code{nmds}, and the last
#'     items dissimilarities for new points. The number of new points
#'     is found as the difference of sizes in \code{nmds} and
#'     \code{dis}.
#' @param neighbors Number of nearest points used to get the starting
#'     locations for new points.
#' @param maxit Maximum number of iterations.
#'
#' @details
#'
#' Function provides an interface to \code{monoMDS} Fortran code with
#' argument \code{NFIX} giving the number of points with fixed
#' coordinates. The number of new points is the difference of number
#' of points in the input ordination (\code{nmds}) and dissimilarities
#' (\code{dis}).
#'
#' @return
#'
#' Function return a list of class \code{"nmds"} (there are no other
#' objects of that type in \pkg{vegan}) with following elements
#' \itemize{
#'   \item{points}{Coordinates of all points, the first being fixed and the
#'     last the new points.}
#'   \item{newpoints}{Coordinates of added new points; the same as the last
#'     ones in the previous item.}
#'   \item{stress}{Stress as estimated here (but needs checking).}
#'   \iters{iters}{Number of iterations.}
#'   \item{cause}{Cause of termination of iterations (numeric).}
#'   \item{seeds}{Starting coordinates for new points.}
#' }
#'
## Open issues:
##
## - Tested only for metaMDS, should be tested also with monoMDS to see
##   if it works there as well.
## - Should accept dis as a symmetric square matrix, and, ultimately, a
##   rectangular matrix of dissimilarities between new and fixed points.
## - Function needs different name!

#' @export
`MDSaddpoints` <-
    function (nmds, dis, neighbors=5, maxit=200)
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


    ## set up initial coordinates as weighted average of nearest
    ## neighbours to old points using 1-diss as weights

    tmp <- matrix(0, newn, ndim)
    for (i in 1:newn) {
        pnt <- order(dis[i,])[seq_len(neighbors)]
        weight <- 1-dis[i,pnt]
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
    isform <- nmds$isform
    ities <- nmds$ities
    iregn <- nmds$iregn
    iscal <- 0L # with NFIX iscal should be 0
    sratmx <- nmds$sratmx
    strmin <- nmds$strmin
    sfgrmn <- nmds$sfgrmn
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
    dim(out$points) <- c(totn, ndim)
    out$newpoints <- out$points[(oldn+1):totn,, drop=FALSE]
    out$seeds <- tmp
    dimnames(out$points)[[1]] <- attr(dis,'Labels')
    dimnames(out$newpoints)[[1]] <- attr(dis,'Labels')[(oldn+1):totn]
    class(out) <- 'nmds'
    out
}
