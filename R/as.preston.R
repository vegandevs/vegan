"as.preston" <-
    function (x, ...) 
{
    if (inherits(x, "preston")) 
        return(x)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    freq <- x[x > 0]
    freq <- ceiling(log2(freq))
    freq <- table(freq)
    nm <- names(freq)
    freq <- as.vector(freq)
    names(freq) <- nm
    class(freq) <- "preston"
    freq
}
