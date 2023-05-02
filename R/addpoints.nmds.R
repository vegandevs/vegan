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
#' @param maxits Maximum number of iterations.
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
`addpoints.nmds` <-
    function (nmds, dis, neighbors=5, maxits=200)
{
## bring list component to local matrix

    points <- nmds$points

## get sizes

    oldn <- nrow(points)
    ndim <- ncol(points)
    totn <- attr(dis,'Size')
    newn <- totn - oldn

## test correspondence

    if (!identical(dimnames(points)[[1]],attr(dis,'Labels')[1:oldn]))
         stop('ordination and dissimilarity matrix do not match')

## decompose disimilarity object to isolate new values

    diss <- as.matrix(dis)[1:oldn,(oldn+1):totn, drop = FALSE]

## set up initial conditions

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

## set up indices

    iidx <- rep((1:oldn),newn)
    jidx <- NULL
    for (i in (oldn+1):totn) jidx <- c(jidx,rep(i,oldn))

## set up ordination parameters.
    nfix <- oldn
    ngrp <-
    istart <- 1
    isform <- nmds$isform
    ities <- nmds$ities
    iregn <- 1
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
                newpoints=matrix(out$points,ncol=ndim)[(oldn+1):totn,, drop=FALSE],
                stress=out$stress,iters=out$iters,cause=out$cause,
                seeds=tmp)
    dimnames(res$points)[[1]] <- attr(dis,'Labels')
    dimnames(res$newpoints)[[1]] <- attr(dis,'Labels')[(oldn+1):totn]
    class(res) <- 'nmds'
    res
}
