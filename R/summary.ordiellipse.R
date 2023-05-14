### Centres and areas of plotted ellipses. The principal axes of the
### conic (oblique ellipse) are found from the eigenvalues of the
### covariance matrix.
`summary.ordiellipse` <-
    function(object, ...)
{
    cnts <- sapply(object, function(x) x$center)
    ## 2nd eigenvalue should be zero if points are on line (like two
    ## points), but sometimes it comes out negative, and area is NaN
    areas <- sapply(object,
                    function(x)
                        sqrt(pmax.int(0, det(x$cov))) * pi * x$scale^2)
    rbind(cnts, `Area` = areas)
}
