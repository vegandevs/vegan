`scores.metaMDS` <-
    function(x, display = c("sites", "species"), shrink = FALSE, choices, ...)
{
    display <- match.arg(display, c("sites","species"), several.ok = TRUE)
    if (missing(choices))
        choices <- seq_len(x$ndim)
    else
        choices <- choices[choices <= x$ndim]
    sites <- NULL
    species <- NULL
    if ("sites" %in% display) {
        sites <- x$points[, choices, drop=FALSE]
        colnames(sites) <- paste0("NMDS", choices)
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
    }
    out <- list()
    out$sites <- sites
    out$species <- species
    if (length(out) == 1)
        out <- out[[1]]
    out
}
