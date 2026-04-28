## Query labels

`labels.envfit` <-
    function(object, ...)
{
    out <- list("vectors" = rownames(object$vectors$arrows),
                "factors" = rownames(object$factors$centroids))
    if (is.null(out$vectors) || is.null(out$factors))
        out <- unlist(out, use.names = FALSE)
    out
}

## Change labels
## THIS WON'T WORK: labels<- IS NOT A GENERIC FUNCTION IN R

`labels<-.envfit` <-
    function(object, value)
{
    if (is.list(value)) { # list of vectors & factors
        if(!is.null(object$vectors) && !is.null(value$vectors))
            rownames(object$vectors$arrows) <- value$vectors
        if (!is.null(object$factors) && !is.null(value$factors))
            rownames(object$factors$centroids) <- value$factors
    } else {
        if(!is.null(object$factors) && !is.null(object$vectors))
            stop("needs a list with both 'vectors' and 'factors' labels")
        if (!is.null(object$factors))
            rownames(object$factors$centroids) <- value
        else
            rownames(object$vectors$arrows) <- value
    }
    object
}
