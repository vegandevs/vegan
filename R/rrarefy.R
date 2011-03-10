### Random rarefied subsample: sample without replacement

`rrarefy` <-
    function(x, sample)
{
    x <- as.matrix(x)
    if (ncol(x) == 1)
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x))
        stop(gettextf(
             "length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length=nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    for (i in 1:nrow(x)) {
        row <- sample(rep(nm, times=x[i,]), sample[i])
        row <- table(row)
        ind <- names(row)
        x[i,] <- 0
        x[i,ind] <- row
    }
    x
}

### Probabilities that species occur in a rarefied 'sample'

`drarefy` <-
    function(x, sample)
{
    if (length(sample) > 1 &&  length(sample) != nrow(x))
        stop(gettextf(
             "length of  'sample' and number of rows of 'x' do not match"))
    x <- drop(as.matrix(x))
    ## dfun is kluge: first item of  vector x must be the sample size,
    ## and the rest  is the community data. This  seemed an easy trick
    ## to evaluate dfun in an apply() instead of a loop.
    dfun <- function(x, sample) {
        J <- sum(x[-1])
        sample <- min(x[1], J)
        1 - exp(lchoose(J - x[-1], sample) - lchoose(J, sample))
    }
    if (length(dim(x)) > 1)
        t(apply(cbind(sample, x), 1, dfun))
    else
        dfun(c(sample, x))
}
