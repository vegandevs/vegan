### Centres and areas of convex hulls (simple polygons).
`summary.ordihull` <-
    function(object, ...)
{
    polyarea <- function(x) {
        n <- nrow(x)
        if (n < 4)
            return(0)
        else
            abs(sum(x[-n,1]*x[-1,2] - x[-1,1]*x[-n,2]))/2
    }
    polycentre <- function(x) {
        n <- nrow(x)
        if (n < 4)
            return(colMeans(x[-n,, drop = FALSE]))
        xy <- x[-n,1]*x[-1,2] - x[-1,1]*x[-n,2]
        A <- sum(xy)/2
        xc <- sum((x[-n,1] + x[-1,1]) * xy)/A/6
        yc <- sum((x[-n,2] + x[-1,2]) * xy)/A/6
        structure(c(xc, yc), names = colnames(x))
    }
    areas <- sapply(object, function(x) polyarea(x))
    cnts <- sapply(object, function(x) polycentre(x))
    rbind(cnts, `Area` = areas)
}
