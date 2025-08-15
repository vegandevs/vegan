`scores.default` <-
    function (x, choices, display = c("sites", "species", "both"),
              tidy = FALSE, ...)
{
    if (is.list(x)) {
        ## why not display <- match.arg(display, names(x)) ?
        display <- match.arg(display)
        if (tidy)
            display <- "both"
        att <- names(x)
    }
    X <- Y <- NULL
    if (is.data.frame(x) && all(sapply(x, is.numeric)))
        x <- as.matrix(x, rownames.force = TRUE)
    if (is.list(x) && display %in% c("sites", "both")) {
        if ("points" %in% att)
            X <- x$points
        else if ("rproj" %in% att)
            X <- x$rproj
        else if ("x" %in% att)
            X <- x$x
        else if ("scores" %in% att)
            X <- x$scores
        else if ("sites" %in% att)
            X <- x$sites
        else if("li" %in% att)
            X <- x$li
        else if("l1" %in% att)
            X <- x$l1
        else stop("cannot find scores")
    }
    if (is.list(x) && display %in% c("species", "both")) {
        if ("species" %in% att)
            Y <- x$species
        else if ("cproj" %in% att)
            Y <- x$cproj
        else if ("rotation" %in% att)
            Y <- x$rotation
        else if ("loadings" %in% att)
            Y <- x$loadings
        else if ("co" %in% att)
            Y <- x$co
        else if ("c1" %in% att)
            Y <- x$c1
        else if (display == "species") # fail if species explicitly requested
            stop("cannot find species scores")
        else { # "both" may be non-chalant: only warn
            warning("cannot find species scores")
        }
    }
    else if (is.numeric(x)) {
        X <- as.matrix(x, rownames.force = TRUE)
        ## as.matrix() changes a score vector to 1-col matrix: this is
        ## a hack which may fail sometimes (but probably less often
        ## than without this hack):

        ## Removed this hack after an issue raised by
        ## vanderleidebastiani in github. He was worried for getting
        ## an error when 'choices' were not given with genuinely 1-dim
        ## (1-col) results. At a second look, it seems that this hack
        ## will fail both with missing 'choices', and also often with
        ## 'choices' given because 'choices' are only applied later,
        ## so that nrow(X) > length(choices). Only vectors (dim arg
        ## missing) should fail here. Let's see...

        ##if (ncol(X) == 1 && nrow(X) == length(choices))
        ##    X <- t(X)
    }
    if (!is.null(X) && NROW(X) && is.null(rownames(X))) {
        rownames(X) <- paste0("site", 1:nrow(X))
    }
    if (!is.null(Y) && NROW(Y) && is.null(rownames(Y))) {
        rownames(Y) <- paste0("spec", 1:nrow(Y))
    }
    if (!is.null(X) && NCOL(X) && is.null(colnames(X))) {
        colnames(X) <- paste0("Dim", 1:ncol(X))
    }
    if (!is.null(Y) && NCOL(Y) && is.null(colnames(Y))) {
        colnames(Y) <- paste0("Dim", 1:ncol(Y))
    }
    if (!missing(choices)) {
        if (!is.null(X))
            X <- X[, choices[choices <= NCOL(X)], drop = FALSE]
        if (!is.null(Y))
            Y <- Y[, choices[choices <= NCOL(Y)], drop = FALSE]
    }
    out <- list("sites" = X, "species" = Y)
    if (tidy) {
        score <- sapply(out, NROW)
        out <- data.frame(do.call(rbind, out),
                          "scores" = rep(names(score), score))
        out$label <- rownames(out)
    }
    if (any(drop <- sapply(out, is.null))) {
        out <- out[!drop]
        if (is.list(out) && length(out) == 1)
            out <- out[[1]]
    }
    out
}
