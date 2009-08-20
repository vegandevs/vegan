### Centres and areas of convex hulls. The area of the *convex* polygon
### is found as the *sum* of the areas of triangles. 
`summary.ordihull` <-
    function(object, ...)
{
    ## The area of triangle from vertices using eq. 8.10 of Spiegel,
    ## Liu & Lipschitz (1999), Mathematical Handbook of Formulas and
    ## Tables (2nd ed.), McGraw & Hill. The hull is closed, so that
    ## first and last point are identical.
    triarea <- function(x) {
        ones <- rep(1, 3)
        if (nrow(x) < 4)
            return(0)
        if (nrow(x) == 4)
            return(abs(det(cbind(x[-1,], ones))/2))
        else {
            sol <- 0
            cnt <- colMeans(x[-1,])
            for (i in 2:nrow(x)) {
                mat <- cbind(rbind(cnt, x[(i-1):i,]), ones)
                sol <- sol + abs(det(mat)/2)
            }
            return(sol)
        }
    }
    areas <- sapply(object, function(x) triarea(x))
    cnts <- sapply(object, function(x) colMeans(x[-1,]))
    rbind(cnts, `Area` = areas)
}
