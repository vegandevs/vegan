`print.poolaccum` <-
    function(x, ...)
{
    rownames(x$means) <- rep("", nrow(x$means))
    print(x$means, ...)
    invisible(x)
}
