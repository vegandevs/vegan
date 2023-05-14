`predict.specaccum` <-
    function(object, newdata, interpolation = c("linear", "spline"), ...)
{
    if (missing(newdata))
        out <- object$richness
    else {
        interpolation <- match.arg(interpolation)
        newdata <- drop(as.matrix(newdata))
        if (length(dim(newdata)) > 1)
            stop("function accepts only one variable as 'newdata'")
        ## Estimation uses lchoose(), but for predict we need to
        ## estimates on non-integer sample sizes and therefore we use
        ## lgamma(). Original "rarefaction" used sample sizes rounded
        ## to integers, but here we can use non-integer data and hence
        ## get different results.
        if (object$method %in% c("exact", "rarefaction")) {
            lg <- function(n, k) {
                ifelse(k <= n, lgamma(pmax.int(n, 0) + 1) - lgamma(k+1) -
                    lgamma(pmax.int(n-k, 0) + 1), -Inf)
            }
            if (object$method == "exact")
                n <- length(object$sites)
            else {
                n <- sum(object$freq)
                newdata <- newdata / length(object$sites) * n
            }
            ldiv <- lg(n, newdata)
            out <- numeric(length(ldiv))
            for(i in seq_along(newdata)) {
                out[i] <- sum(1 - exp(lg(n-object$freq, newdata[i])
                                      - ldiv[i]))
            }
        } else if (object$method == "coleman") {
            ## "coleman" also works on non-integer newdata
            n <- length(object$sites)
            out <- sapply(newdata,
                          function(x) sum(1 - (1 - x/n)^object$freq))
        } else {
            ## Other methods do not accept non-integer newdata, but we
            ## can interpolate
            if (interpolation == "linear")
                out <- approx(x = object$sites, y = object$richness,
                              xout = newdata, rule = 1)$y
            else
                out <- spline(x = object$sites, y = object$richness,
                              xout = newdata, ...)$y
        }
    }
    out
}
