`scores.envfit` <-
    function (x, display, choices, ...)
{
    display <- match.arg(display,
                         c("vectors", "bp", "factors", "cn"),
                         several.ok = TRUE)
    out <- list()
    if (any(display %in% c("vectors", "bp"))) {
        vects <- x$vectors$arrows[, , drop = FALSE]
        if (!missing(choices))
            vects <- vects[, choices, drop=FALSE]
        if (!is.null(vects))
            out$vectors <- sqrt(x$vectors$r) * vects
    }
    if (any(display %in% c("factors", "cn"))) {
        facts <- x$factors$centroids[, , drop = FALSE]
        if (!missing(choices))
            facts <- facts[, choices, drop=FALSE]
        out$factors <- facts
    }
    if (length(out) == 1)
        out <- out[[1]]
    out
}

