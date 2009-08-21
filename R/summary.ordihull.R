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
    areas <- sapply(object, function(x) polyarea(x))
    cnts <- sapply(object, function(x) colMeans(x[-1,]))
    rbind(cnts, `Area` = areas)
}
