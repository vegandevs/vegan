### Function to find distances of points to centroids

`centredist` <- function(x) UseMethod("centredist")

### Function to find distances from each sampling unit to each centroid
### for vegan::betadisper result.
### x (input): result object from vegan::betadisper
###
### Originally published in
### https://stackoverflow.com/questions/77391007/ and in github issue
### #606

`centredist.betadisper` <-
    function(x)
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
    class(out) <- "centredist"
    out
}

`print.centredist`<-
    function(x, ...)
{
     print(data.frame("centre" = x$centre, "nearest" =  x$nearest, x$distances))
     invisible(x)
}

