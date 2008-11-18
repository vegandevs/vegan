`print.capscale` <-
    function(x, ...)
{
    NextMethod("print", x, ...)
    if (!is.null(x$metaMDSdist))
        cat("metaMDSdist transformed data:", x$metaMDSdist, "\n\n") 
}
