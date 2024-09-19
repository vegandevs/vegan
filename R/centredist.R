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
    function(x, centres, distance = c("euclidean", "mahalanobis"),
             display = "sites", w, ...)
{
    centres <- factor(centres)
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
        if (min(table(centres)) <= ncol(x))
        if (max(table(centres)) <= ncol(x))
            stop(gettextf(
                "no group is larger than the number of dimensions (%d)",
                ncol(x)))
        if (min(table(centres)) <= ncol(x))
            warning(
                gettextf("groups smaller or equal to no. of dimensions (%d) will be NA",
                         ncol(x)))
        for(cl in levels(centres)) {
            cv[[cl]] <- cov.wt(x[centres == cl,, drop=FALSE],
                               wt = w[centres == cl])
        }
        dis <-
            sapply(levels(centres),
                   function (cl) if (det(cv[[cl]]$cov) <= .Machine$double.eps)
                                     rep(NA, nrow(x))
                                 else
                                     mahalanobis(x, cv[[cl]]$center,
                                                 cv[[cl]]$cov))
    } else { # euclidean
        sumw <- sum(w)
        w <- w/sumw
        cnt <- sapply(levels(centres), function(cl)
            apply(x[centres == cl,, drop=FALSE], 2,
                  weighted.mean, w = w[centres == cl]))
        dis <- apply(cnt, 2, function(z) w * rowSums(sweep(x, 2, z)^2))
        dis <- sqrt(sumw * dis)
    }
    nearest <- levels(centres)[apply(dis, 1, which.min)]
    out <- list("centre" = centres, "nearest" = nearest, "distances" = dis)
    attr(out, "distance") <- distance
    class(out) <- "centredist"
    out
}

## Above: dist="eucl": sum of squared distances to their centre should
## give eigenvalues for consistent estimates. Consistent distances to
## rda model 'm': using scaling="site" with const =
## sqrt(m$tot.chi*nobs(m)) and w = rep(1/nobs(m), nobs(m)). Consistent
## to wchmdscale with const = sqrt(m$tot.chi*nobs(m)) and w = rep(1,
## nobs(m)). The cca results can be used directly (naturally, with
## scaling="site").

## For the Americans

`centerdist` <- function(...) centredist(...)

`print.centredist`<-
    function(x, ...)
{
    cat("\ndistances: ", attr(x, "distance"), "\n\n")
    print(data.frame("centre" = x$centre, "nearest" =  x$nearest, x$distances))
    invisible(x)
}

