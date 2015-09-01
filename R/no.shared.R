"no.shared" <-
    function(x)
{
    N <- nrow(x <- as.matrix(x))
    d <- .C("veg_distance", x = as.double(x), nr = N, nc = ncol(x),
            d = double(N * (N - 1)/2), diag = as.integer(FALSE),
            method = as.integer(99), PACKAGE="vegan")$d
    d <- as.logical(d)
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "method") <- "no.shared"
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    d        
}
