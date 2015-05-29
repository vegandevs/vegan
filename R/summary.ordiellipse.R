### Centres and areas of plotted ellipses. The principal axes of the
### conic (oblique ellipse) are found from the eigenvalues of the
### covariance matrix.
`summary.ordiellipse` <-
    function(object, ...)
{
    cnts <- sapply(object, function(x) x$center)
    ## two points are exactly on line and have zero area ellipses
    areas <- sapply(object,
                    function(x) if (x$n.obs <= 2) 0 else
                    prod(sqrt(eigen(x$cov)$values)) * pi * x$scale^2)
    rbind(cnts, `Area` = areas)
}
