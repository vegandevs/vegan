`as.fisher` <-
    function (x, ...) 
{
    if (inherits(x, "fisher")) 
        return(x)
    ## is not fisher but a 1 x n data.frame or matrix: matrix is faster
    x <- as.matrix(x)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    freq <- x[x > 0]
    freq <- table(freq, deparse.level = 0)
    nm <- names(freq)
    freq <- as.vector(freq)
    names(freq) <- nm
    class(freq) <- "fisher"
    freq
}
