`veganCovEllipse` <-
    function(cov, center = c(0,0), scale = 1, npoints = 100)
{
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    t(center + scale * t(Circle %*% chol(cov)))
}
