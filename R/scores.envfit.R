`scores.envfit` <-
    function (x, display, choices, ...)
{
    display <- match.arg(display, c("vectors", "bp", "factors",
        "cn"))
    if (display %in% c("vectors", "bp")) {
        out <- x$vectors$arrows[, , drop = FALSE]
        if (!is.null(out))
            out <- sweep(out, 1, sqrt(x$vectors$r), "*")
    }
    else out <- x$factors$centroids[, , drop = FALSE]
    if (!missing(choices) && !is.null(out)) {
        if (length(choices) <= ncol(out))
            out <- out[, choices, drop = FALSE]
        else
            stop(gettextf(
                "you requested  %d dimensions, but 'envfit' only has %d",
                length(choices), ncol(out)))
    }
    out
}

