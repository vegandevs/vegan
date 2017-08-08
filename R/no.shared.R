`no.shared` <-
    function(x)
{
    x <- as.matrix(x)
    d <- .Call(do_vegdist, x, as.integer(99))
    d <- as.logical(d)
    attr(d, "Size") <- NROW(x)
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "method") <- "no.shared"
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    d
}
