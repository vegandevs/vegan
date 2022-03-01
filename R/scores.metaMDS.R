`scores.metaMDS` <-
    function(x, display = c("sites", "species"), shrink = FALSE, choices,
             tidy = FALSE, ...)
{
    display <- match.arg(display, c("sites","species"), several.ok = TRUE)
    if (missing(choices))
        choices <- seq_len(x$ndim)
    else
        choices <- choices[choices <= x$ndim]
    out <- list()
    if ("sites" %in% display) {
        sites <- x$points[, choices, drop=FALSE]
        colnames(sites) <- paste0("NMDS", choices)
        out$sites <- sites
    }
    if ("species" %in% display && !is.null(x$species)) {
        species <- x$species[, choices, drop=FALSE]
        colnames(species) <- paste0("NMDS", choices)
        if (shrink) {
            mul <- sqrt(attr(X, "shrinkage"))
            if (is.null(mul))
                warning("species cannot be shrunken, because they were not expanded")
            else {
                mul <- mul[choices]
                cnt <- attr(X, "centre")
                X <- sweep(X, 2, cnt, "-")
                X <- sweep(X, 2, mul, "*")
                X <- sweep(X, 2, cnt, "+")
            }
        }
        out$species <- species
    }
    if (tidy) {
        if (length(out) == 0) # no scores (species scores may not exist)
            return(NULL)
        group <- sapply(out, nrow)
        group <- rep(names(group), group)
        out <- do.call(rbind, out)
        label <- rownames(out)
        out <- as.data.frame(out)
        out$score <- group
        out$label <- label
    }
    ## only two kind of scores, return NULL, matrix, or a list of scores
    switch(length(out), out[[1]], out)
}
