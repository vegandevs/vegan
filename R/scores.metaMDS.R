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
    if ("species" %in% display && !is.null(x$species) && !all(is.na(x$species))) {
        species <- x$species[, choices, drop=FALSE]
        colnames(species) <- paste0("NMDS", choices)
        if (shrink) {
            ## [,choices] drops attributes
            mul <- sqrt(attr(x$species, "shrinkage"))
            cnt <- attr(x$species, "centre")
            if (is.null(mul))
                message("species are not shrunken, because they were not expanded")
            else {
                mul <- mul[choices]
                cnt <- cnt[choices]
                species <- sweep(species, 2, cnt, "-")
                species <- sweep(species, 2, mul, "*")
                species <- sweep(species, 2, cnt, "+")
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
    if (length(out) == 1)
        out[[1]]
    else
        out
}
