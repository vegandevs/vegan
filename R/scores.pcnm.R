`scores.pcnm` <-
    function(x, choices, ...)
{
    if (missing(choices))
        x$vectors
    else
        x$vectors[, choices]
}
