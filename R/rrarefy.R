### Random rarefied subsample: sample without replacement

`rrarefy` <-
    function(x, sample, iterations = 1)
{
    if (!identical(all.equal(x, round(x)), TRUE)) 
        stop("function is meaningful only for integers (counts)")
    x <- as.matrix(x)
    if (ncol(x) == 1) 
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x)) 
        stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length = nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    if (any(rowSums(x) < sample)) 
        warning("Some row sums < 'sample' and are not rarefied")
    for (i in 1:nrow(x)) {
        if (sum(x[i, ]) <= sample[i]) 
            next
        row <- sample(rep(nm, times = x[i, ]), sample[i])
        row <- table(row)
        ind <- names(row)
        y <- x
        y[i, ] <- 0
        y[i, ind] <- row
        counter <- 1
        if(iterations > 1) {
            z <- x
            itr <- iterations-1
            for(k in c(1:itr)) {
                row <- sample(rep(nm, times = x[i, ]), sample[i])
                row <- table(row)
                ind <- names(row)
                z[i, ] <- 0
                z[i, ind] <- row
                z[i,] <- y[i,] + (1/(counter + 1)) * (z[i,] - y[i,])
                y[i,] <- z[i,]
                counter <- counter + 1
            }
        }
    }
    y
}

### Probabilities that species occur in a rarefied 'sample'

`drarefy` <-
    function(x, sample)
{
    if (!identical(all.equal(x, round(x)), TRUE)) 
        stop("function accepts only integers (counts)")
    if (length(sample) > 1 &&  length(sample) != nrow(x))
        stop(gettextf(
             "length of 'sample' and number of rows of 'x' do not match"))
    x <- drop(as.matrix(x))
    ## warn on too large samples
    if (is.matrix(x))
        rs <- rowSums(x)
    else
        rs <- sum(x)
    if (any(rs) < sample)
        warning("Some row sums < 'sample' and probabilities either 0 or 1")
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
