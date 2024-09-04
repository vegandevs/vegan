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
    function(x, centres, distance = c("euclidean", "mahalanobis"), ...)
{
    ## assume so far that x is a matrix that can be handled as such
    centres <- factor(centres)
    distance <- match.arg(distance)
    cv <- list()
    if (distance == "mahalanobis") {
        for(cl in levels(centres)) {
            cv[[cl]] <- cov.wt(x[centres == cl,, drop=FALSE])
        }
        dis <-
            sapply(levels(centres),
                   function (cl) if (det(cv[[cl]]$cov) <= .Machine$double.eps)
                                     rep(NA, nrow(x))
                                 else
                                     mahalanobis(x, cv[[cl]]$center,
                                                 cv[[cl]]$cov))
    } else { # euclidean
        cnt <- sapply(levels(centres),
                      function(cl) colMeans(x[centres == cl,, drop=FALSE]))
        dis <- apply(cnt, 2, function(z) rowSums(sweep(x, 2, z)^2))
        dis <- sqrt(dis)
    }
    nearest <- levels(centres)[apply(dis, 1, which.min)]
    out <- list("centre" = centres, "nearest" = nearest, "distances" = dis)
    attr(out, "distance") <- distance
    class(out) <- "centredist"
    out
}

`print.centredist`<-
    function(x, ...)
{
    cat("\ndistances: ", attr(x, "distance"), "\n\n")
    print(data.frame("centre" = x$centre, "nearest" =  x$nearest, x$distances))
    invisible(x)
}

