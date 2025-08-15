### Random rarefied subsample: sample without replacement

`rrarefy` <-
    function(x, sample)
{
    x <- as.matrix(x, rownames.force = TRUE)
    if (!isTRUE(all.equal(x, round(x))))
        stop("function is meaningful only for integers (counts)")
    ## x may not be exactly integer, since, e.g., sqrt(2)^2 != 2
    if (!is.integer(x))
        x <- round(x)
    minobs <- min(x[x > 0])
    if (minobs > 1)
        warning(gettextf("function should be used for observed counts, but smallest count is %d", minobs))
    if (ncol(x) == 1)
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x))
        stop(gettextf(
             "length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length=nrow(x))
    ## warn if something cannot be rarefied
    if (any(rowSums(x) < sample))
        warning("some row sums < 'sample' and are not rarefied")
    for (i in 1:nrow(x)) {
        x[i,] <- .Call(do_rrarefy, x[i,], sample[i], PACKAGE = "vegan")
    }
    x
}

### Probabilities that species occur in a rarefied 'sample'

`drarefy` <-
    function(x, sample)
{
    if (!isTRUE(all.equal(x, round(x))))
        stop("function accepts only integers (counts)")
    x <- round(x)
    minobs <- min(x[x > 0])
    if (minobs > 1)
        warning(gettextf("most observed count data have counts 1, but smallest count is %d", minobs))
    if (length(sample) > 1 &&  length(sample) != nrow(x))
        stop(gettextf(
             "length of 'sample' and number of rows of 'x' do not match"))
    x <- drop(as.matrix(x))
    ## warn on too large samples
    if (is.matrix(x))
        rs <- rowSums(x)
    else
        rs <- sum(x)
    if (any(rs < sample))
        warning("some row sums < 'sample' and probabilities either 0 or 1")
    ## dfun is kluge: first item of  vector x must be the sample size,
    ## and the rest  is the community data. This  seemed an easy trick
    ## to evaluate dfun in an apply() instead of a loop.
    dfun <- function(x) {
        J <- sum(x[-1])
        sample <- min(x[1], J)
        1 - exp(lchoose(J - x[-1], sample) - lchoose(J, sample))
    }
    if (length(dim(x)) > 1)
        t(apply(cbind(sample, x), 1, dfun))
    else
        dfun(c(sample, x))
}
