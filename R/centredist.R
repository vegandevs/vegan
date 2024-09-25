### Function to find distances of points to centroids

`centredist` <- function(x, ...) UseMethod("centredist")

### Function to find distances from each sampling unit to each centroid
### for vegan::betadisper result.
### x (input): result object from vegan::betadisper
###
### Originally published in
### https://stackoverflow.com/questions/77391007/ and in github issue
### #606

`centredist.betadisper` <-
    function(x, ...)
{
    cnt <- x$centroids
    coord <- x$vectors
    pos <- which(x$eig >= 0)
    neg <- which(x$eig < 0)
    d <- apply(cnt[,pos], 1,
               function(z) rowSums(sweep(coord[,pos], 2, z)^2))
    if (length(neg))
        d <- d - apply(cnt[, neg], 1,
                       function(z) rowSums(sweep(coord[,neg], 2, z)^2))
    d <- as.data.frame(sqrt(d))
    nearest <- levels(x$group)[apply(d, 1, which.min)]
    out <- list("centre" = x$group, "nearest" = nearest, "distances" = d)
    attr(out, "distance") <- attr(x, "method")
    class(out) <- "centredist"
    out
}

`centredist.default` <-
    function(x, group, distance = c("euclidean", "mahalanobis"),
             display = "sites", w, ...)
{
    group <- factor(group)
    distance <- match.arg(distance)
    if (missing(w))
        w <- if (is.atomic(x)) attr(x, "weights")
             else weights(x, display = display, ...)
    x <- scores(x, display = display, ...)
    if (is.null(w))
        w <- rep(1, nrow(x))
    if (inherits(x, c("dbrda", "wcmdscale")) && any(eigenvals(x) < 0))
        warning("axes with negative eigenvalues are ignored")
    cv <- list()
    if (distance == "mahalanobis") {
        if (min(table(group)) <= ncol(x))
        if (max(table(group)) <= ncol(x))
            stop(gettextf(
                "no group is larger than the number of dimensions (%d)",
                ncol(x)))
        if (min(table(group)) <= ncol(x))
            warning(
                gettextf("groups smaller or equal to no. of dimensions (%d) will be NA",
                         ncol(x)))
        for(cl in levels(group)) {
            cv[[cl]] <- cov.wt(x[group == cl,, drop=FALSE],
                               wt = w[group == cl])
        }
        dis <-
            sapply(levels(group),
                   function (cl) if (det(cv[[cl]]$cov) <= .Machine$double.eps)
                                     rep(NA, nrow(x))
                                 else
                                     mahalanobis(x, cv[[cl]]$center,
                                                 cv[[cl]]$cov))
    } else { # euclidean
        sumw <- sum(w)
        w <- w/sumw
        cnt <- sapply(levels(group), function(cl)
            apply(x[group == cl,, drop=FALSE], 2,
                  weighted.mean, w = w[group == cl]))
        dis <- apply(cnt, 2, function(z) w * rowSums(sweep(x, 2, z)^2))
        dis <- sqrt(sumw * dis)
    }
    nearest <- levels(group)[apply(dis, 1, which.min)]
    out <- list("centre" = group, "nearest" = nearest, "distances" = dis)
    attr(out, "distance") <- distance
    class(out) <- "centredist"
    out
}

## Above: dist="eucl": sum of squared distances to their centre should
## give eigenvalues for consistent estimates. For consistent
## distances, see fucntions centredist.cca and centredist.rda below.


## cca: use scaling by scores type
`centredist.cca` <-
    function(x, group, distance = c("euclidean", "mahalanobis"),
             display = c("sites", "species"), rank = 2, ...)
{
    ## only accept displays "sites", "species"
    display <- match.arg(display)
    if (pmatch(rank, "full", nomatch = FALSE))
        rank <- max(x$CCA$rank, 0) + max(x$CA$rank, 0)
    choices = seq_len(rank)
    centredist.default(x = x, group = group, distance = distance,
                       display = display, scaling = display, choices = choices,
                       ...)
}

## rda, dbrda, capscale: scaling by scores type, const and w give
## results consistent with the returned eigenvalue
`centredist.rda` <-
        function(x, group, distance = c("euclidean", "mahalanobis"),
                 display = "sites", rank = 2, ...)
{
    display <- match.arg(display)
    if (pmatch(rank, "full", nomatch = FALSE))
        rank <- max(x$CCA$rank, 0) + max(x$CA$rank, 0)
    choices = seq_len(rank)
    N <- nobs(x)
    w <- weights(x, displaya = display)
    centredist.default(x = x, group = group, distance = distance,
                       display = display, scaling = display,
                       choices = choices,
                       const = sqrt(x$tot.chi * N),
                       w = w/sum(w), ...)
}

## wcmdscale class (i.e., with eig=TRUE) may have negative eigenvalues

`centredist.wcmdscale` <-
    function(x, group, distance = c("euclidean", "mahalanobis"),
             display = "sites", rank = 2, ...)
{
    display <- match.arg(display)
    distance <- match.arg(distance)
    if (is.numeric(rank)) {
        centredist.default(x = x, group = group, distance = distance,
                           display = display, choices = seq_len(rank),
                           ...)
    } else if (pmatch(rank, "full", nomatch = FALSE)) {
        if (distance == "mahalanobis")
            stop("rank='full' is only possible with distance='euclidean'")
        xre <- centredist(x = x$points, group = group, distance = distance,
                          ...)
        if (!is.null(x$negaxes)) {
            xim <- centredist(x = x$negaxes, group = group, distance = distance,
                              ...)$distances
            xre$distances <- sqrt(xre$distances^2 - xim^2)
            xre
        }
    } else (stop("'rank' should be integer of 'full'"))
}

## For the Americans

`centerdist` <- function(...) centredist(...)

`print.centredist`<-
    function(x, ...)
{
    cat("\ndistances: ", attr(x, "distance"), "\n\n")
    print(data.frame("centre" = x$centre, "nearest" =  x$nearest, x$distances))
    invisible(x)
}

