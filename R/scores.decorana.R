`scores.decorana` <-
    function (x, display="sites", choices = 1:4, origin=TRUE,
              tidy = FALSE, ...)
{
    display <- match.arg(display, c("sites", "species", "both"), several.ok = TRUE)
    ## return "both" in tidy scores
    if(tidy)
        display <- "both"
    out <- list()
    if(any(c("sites", "both") %in% display)) {
        sites <- x$rproj
        if (origin)
            sites <- sweep(sites, 2, x$origin, "-")
        out$sites <- sites[, choices, drop=FALSE]
    }
    if(any(c("species", "both") %in% display)) {
        species <- x$cproj
        if (origin)
            species <- sweep(species, 2, x$origin, "-")
        out$species <- species[, choices]
    }
    if (tidy) {
        if (length(out) == 0) # no scores (never TRUE?)
            return(NULL)
        group <- sapply(out, nrow)
        group <- rep(names(group), group)
        out <- do.call(rbind, out)
        label <- rownames(out)
        out <- as.data.frame(out)
        out$score <- group
        out$label <- label
        wts <- rep(NA, nrow(out))
        if (any(take <- group == "sites"))
            wts[take] <- weights(x, display="sites")
        if (any(take <- group == "species"))
            wts[take] <- weights(x, display="species")
        out$weight <- wts
    }
    ## two kind of scores: return NULL, matrix or a list
    if (length(out) == 1)
        out[[1]]
    else
        out
}
