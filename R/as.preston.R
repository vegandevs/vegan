`as.preston` <-
    function (x, tiesplit = TRUE, ...)
{
    if (inherits(x, "preston"))
        return(x)
    ## practically integer
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    ## need exact integers, since, e.g., sqrt(2)^2 - 2 = 4.4e-16 and
    ## tie breaks fail
    if (!is.integer(x))
        x <- round(x)
    x <- x[x > 0]
    if (tiesplit) {
        ## Assume log2(2^k) == k *exactly* for integer k
        xlog2 <- log2(x)
        ties <- xlog2 == ceiling(xlog2)
        tiefreq <- table(xlog2[ties])
        notiefreq <- table(ceiling(xlog2[!ties]))
        itie <- as.numeric(names(tiefreq)) + 1
        nitie <- as.numeric(names(notiefreq)) + 1
        freq <- numeric(max(itie+1, nitie))
        ## split tied values between two adjacent octaves
        freq[itie] <- tiefreq/2
        freq[itie+1] <- freq[itie+1] + tiefreq/2
        freq[nitie] <- freq[nitie] + notiefreq
    } else {
        xlog2 <- ceiling(log2(x))
        tmp <- table(xlog2)
        indx <- as.numeric(names(tmp)) + 1
        freq <- numeric(max(indx))
        freq[indx] <- tmp
    }
    names(freq) <- seq_along(freq) - 1
    ## remove empty octaves
    freq <- freq[freq>0]
    class(freq) <- "preston"
    freq
}
