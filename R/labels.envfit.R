## Query labels used in plotting

`labels.envfit` <-
    function(object, ...)
{
    out <- list("vectors" = rownames(object$vectors$arrows),
                "factors" = rownames(object$factors$centroids))
    if (is.null(out$vectors) || is.null(out$factors))
        out <- unlist(out, use.names = FALSE)
    out
}

## Change names of variables in the envfit result

`names<-.envfit` <-
    function(x, value)
{
    if (is.list(value)) { # list of vectors & factors
        if(!is.null(x$vectors) && !is.null(value$vectors))
            rownames(x$vectors$arrows) <- value$vectors
        if (!is.null(x$factors) && !is.null(value$factors))
            rownames(x$factors$centroids) <- value$factors
    } else {
        if(!is.null(x$factors) && !is.null(x$vectors))
            stop("needs a list with both 'vectors' and 'factors' labels")
        if (!is.null(x$factors))
            rownames(x$factors$centroids) <- value
        else
            rownames(x$vectors$arrows) <- value
    }
    x
}
