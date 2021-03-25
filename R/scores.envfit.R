`scores.envfit` <-
    function (x, display, choices, arrow.mul = 1, tidy = FALSE, ...)
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
            out$vectors <- arrow.mul * sqrt(x$vectors$r) * vects
    }
    if (any(display %in% c("factors", "cn"))) {
        facts <- x$factors$centroids[, , drop = FALSE]
        if (!missing(choices))
            facts <- facts[, choices, drop=FALSE]
        out$factors <- facts
    }
    if (tidy) {
        if (length(out) == 0) # no scores
            return(NULL)
        group <- sapply(out, nrow)
        group <- rep(names(group), group)
        out <- do.call(rbind, out)
        label <- rownames(out)
        out <- as.data.frame(out)
        out$score <- group
        out$label <- label
    }
    ## only two kind of scores: return NULL, matrix or a list
    switch(length(out), out[[1]], out)
}

