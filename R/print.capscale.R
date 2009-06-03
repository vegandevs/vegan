`print.capscale` <-
    function(x, ...)
{
    NextMethod("print", x, ...)
    if (any(x$CA$eig < 0))
        cat("NB. Some eigenvalues are negative\n\n")
    if (!is.null(x$metaMDSdist))
        cat("metaMDSdist transformed data:", x$metaMDSdist, "\n\n") 
}
