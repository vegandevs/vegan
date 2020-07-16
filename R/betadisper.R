`betadisper` <-
    function(d, group, type = c("median","centroid"), bias.adjust=FALSE,
             sqrt.dist = FALSE, add = FALSE)
{
    ## inline function for double centring. We used .C("dblcen", ...,
    ## PACKAGE = "stats") which does not dublicate its argument, but
    ## it was removed from R in r60360 | ripley | 2012-08-22 07:59:00
    ## UTC (Wed, 22 Aug 2012) "more conversion to .Call, clean up".
    dblcen <- function(x, na.rm = TRUE) {
        cnt <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2L, cnt, check.margin = FALSE)
        cnt <- rowMeans(x, na.rm = na.rm)
        sweep(x, 1L, cnt, check.margin = FALSE)
    }
    ## inline function for spatial medians
    spatialMed <- function(vectors, group, pos) {
        axes <- seq_len(NCOL(vectors))
        spMedPos <- ordimedian(vectors, group, choices = axes[pos])
        spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
        cbind(spMedPos, spMedNeg)
    }
    ## inline function for centroids
    centroidFUN <- function(vec, group) {
        cent <- apply(vec, 2,
                      function(x, group) tapply(x, INDEX = group, FUN = mean),
                      group = group)
        if(!is.matrix(cent)) { ## if only 1 group, cent is vector
            cent <- matrix(cent, nrow = 1,
                           dimnames = list(as.character(levels(group)),
                           paste0("Dim", seq_len(NCOL(vec)))))
        }
        cent
    }
    ## inline function for distance computation
    Resids <- function(x, c) {
        if(is.matrix(c))
            d <- x - c
        else
            d <- sweep(x, 2, c)
        rowSums(d^2)
    }
    ## Tolerance for zero Eigenvalues
    TOL <- sqrt(.Machine$double.eps)
    ## uses code from stats:::cmdscale by R Core Development Team
    if(!inherits(d, "dist"))
        stop("distances 'd' must be a 'dist' object")
    ## Someone really tried to analyse correlation like object in range -1..+1
    if (any(d < -TOL, na.rm = TRUE))
        stop("dissimilarities 'd' must be non-negative")
    ## adjust to avoid negative eigenvalues (if they disturb you)
    if (sqrt.dist)
        d <- sqrt(d)
    if (is.logical(add) && isTRUE(add))
        add <- "lingoes"
    if (is.character(add)) {
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- addLingoes(as.matrix(d))
            d <- sqrt(d^2 + 2 * ac)
        }
        else if (add == "cailliez") {
            ac <- addCailliez(as.matrix(d))
            d <- d + ac
        }
    }
    if(missing(type))
        type <- "median"
    type <- match.arg(type)
    ## checks for groups - need to be a factor for later
    group <- if(!is.factor(group)) {
        as.factor(group)
    } else { ## if already a factor, drop empty levels
        droplevels(group, exclude = NA) # need exclude = NA under Rdevel r71113
    }
    n <- attr(d, "Size")
    x <- matrix(0, ncol = n, nrow = n)
    x[row(x) > col(x)] <- d^2
    ## site labels
    labs <- attr(d, "Labels")
    ## remove NAs in group
    if(any(gr.na <- is.na(group))) {
        group <- group[!gr.na]
        x <- x[!gr.na, !gr.na]
        ## update n otherwise C call crashes
        n <- n - sum(gr.na)
        ## update labels
        labs <- labs[!gr.na]
        message("missing observations due to 'group' removed")
    }
    ## remove NA's in d
    if(any(x.na <- apply(x, 1, function(x) any(is.na(x))))) {
        x <- x[!x.na, !x.na]
        group <- group[!x.na]
        ## update n otherwise C call crashes
        n <- n - sum(x.na)
        ## update labels
        labs <- labs[!x.na]
        message("missing observations due to 'd' removed")
    }
    x <- x + t(x)
    x <- dblcen(x)
    e <- eigen(-x/2, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    ## Remove zero eigenvalues
    eig <- eig[(want <- abs(eig) > max(TOL, TOL * eig[1L]))]
    ## scale Eigenvectors
    vectors <- vectors[, want, drop = FALSE] %*% diag(sqrt(abs(eig)),
                               nrow = length(eig))
    ## store which are the positive eigenvalues
    pos <- eig > 0
    ## group centroids in PCoA space
    centroids <-
        switch(type,
               centroid = centroidFUN(vectors, group),
               median = spatialMed(vectors, group, pos)
               )
    ## for each of the groups, calculate distance to centroid for
    ## observation in the group
    ## Uses in-line Resids function as we want LAD residuals for
    ## median method, and LSQ residuals for centroid method
    dist.pos <- Resids(vectors[, pos, drop=FALSE],
                       centroids[group, pos, drop=FALSE])
    dist.neg <- 0
    if(any(!pos))
        dist.neg <- Resids(vectors[, !pos, drop=FALSE],
                           centroids[group, !pos, drop=FALSE])

    ## zij are the distances of each point to its group centroid
    if (any(dist.neg > dist.pos)) {
        ## Negative squared distances give complex valued distances:
        ## take only the real part (which is zero). Github issue #306.
        warning("some squared distances are negative and changed to zero")
        zij <- Re(sqrt(as.complex(dist.pos - dist.neg)))
    } else {
        zij <- sqrt(dist.pos - dist.neg)
    }
    if (bias.adjust) {
        n.group <- as.vector(table(group))
        zij <- zij*sqrt(n.group[group]/(n.group[group]-1))
    }
    ## add in correct labels
    if (any(want))
        colnames(vectors) <- names(eig) <-
            paste("PCoA", seq_along(eig), sep = "")
    if(is.matrix(centroids))
        colnames(centroids) <- names(eig)
    else
        names(centroids) <- names(eig)
    rownames(vectors) <- names(zij) <- labs
    retval <- list(eig = eig, vectors = vectors, distances = zij,
                   group = group, centroids = centroids, call = match.call())
    class(retval) <- "betadisper"
    attr(retval, "method") <- attr(d, "method")
    attr(retval, "type") <- type
    attr(retval, "bias.adjust") <- bias.adjust
    retval
}
