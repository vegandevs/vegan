`labels.envfit` <-
    function(object, ...)
{
    out <- list("vectors" = rownames(object$vectors$arrows),
                "factors" = rownames(object$factors$centroids))
    if (is.null(out$vectors) || is.null(out$factors))
        out <- unlist(out, use.names = FALSE)
    out
}
