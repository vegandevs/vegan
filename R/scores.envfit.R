`scores.envfit` <-
    function (x, display = c("vectors", "factors"), choices,
              arrow.mul = 1, tidy = FALSE, ...)
{
    display <- match.arg(display,
                         c("vectors", "bp", "factors", "cn"),
                         several.ok = TRUE)
    if (tidy)
        display <- c("vectors", "factors")
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
        ## out <- data.frame("score" = factor(group), label=label,
        ##                  out, row.names = NULL)

        ## Command above would be compatible with ordination fortify,
        ## but current ggvegan::fortify.envfit is inconsistent
        out <- data.frame(label,
                          type = factor(group, levels = c("factors", "vectors"),
                                        labels = c("Centroid", "Vector")),
                          out, row.names = NULL)
    }
    ## Return NULL, matrix of one type of scores, list or tidy data frame
    switch(as.character(length(out)),
           "0" = NULL,
           "1" = out[[1]],
           out)
}

