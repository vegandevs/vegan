"scores.metaMDS" <-
    function(x, display = c("sites", "species"), shrink = FALSE, choices, ...)
{
    display <- match.arg(display)
    if (missing(choices))
        choices <- 1:x$dims
    if (display == "sites")
        X <- x$points
    else if (display == "species") {
        X <- x$species
        if (shrink) {
            mul <- sqrt(attr(X, "shrinkage"))
            if (is.null(mul))
                warning("Species cannot be shrunken, because they were not expanded")
            else {
                cnt <- attr(X, "centre")
                X <- sweep(X, 2, cnt, "-")
                X <- sweep(X, 2, mul, "*")
                X <- sweep(X, 2, cnt, "+")
            }
        }
    }
    colnames(X) <- paste("NMDS", 1:ncol(X), sep="")
    X[, choices]
}

