#' Slope of Rarefunction Curve at Given Sample Size
#'
#' Function evaluates the derivative of the rarefaction
#' function at given sample size. The derivative was
#' directly derived from the expression used in \code{rarefy}.
#'
#' @param x Community counts, either an integer vector for a single
#' site or a data frame or matrix with each row giving site vectors.
#' @param sample Sample sizes where the derivatives are evaluated; can
#' be real
#'
`rareslope` <-
    function(x, sample)
{
    ## 'x' must be integers ('sample' need not be)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("community data 'x' must be integers (counts)")
    slope <- function(x, sample) {
        x <- x[x>0]
        J <- sum(x)
        ## Replace Hurlbert's factorials with gamma() functions and do
        ## some algebra for derivatives. NB., rarefy() does not use
        ## factorials but lchoose directly.
        d <- digamma(pmax(J-sample+1, 1)) - digamma(pmax(J-x-sample+1, 1))
        g <- lgamma(pmax(J-x+1, 1)) + lgamma(pmax(J-sample+1, 1)) -
            lgamma(pmax(J-x-sample+1, 1)) - lgamma(J+1)
        d <- d*exp(g)
        sum(d[is.finite(d)])
    }
    if (length(dim(x)) == 2)
        out <- sapply(sample, function(n) apply(x, 1, slope, sample = n))
    else
        out <- sapply(sample, function(n) slope(x, sample=n))
    out <- drop(out)
    if (length(sample) > 1) {
        if (is.matrix(out))
            colnames(out) <- paste0("N", sample)
        else
            names(out) <- paste0("N", sample)
    }
    out
}
