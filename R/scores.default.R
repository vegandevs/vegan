"scores.default" <-
    function (x, choices, display = c("sites", "species"), ...) 
{
    display <- match.arg(display)
    att <- names(x)
    if (is.data.frame(x) && all(sapply(x, is.numeric)))
        x <- as.matrix(x)
    if (is.list(x) && display == "sites") {
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
        else stop("Can't find scores")
    }
    else if (is.list(x) && display == "species") {
        if ("species" %in% att)
            X <- x$species
        else if ("cproj" %in% att) 
            X <- x$cproj
        else if ("rotation" %in% att) 
            X <- x$rotation
        else if ("loadings" %in% att) 
            X <- x$loadings
        else if ("co" %in% att)
            X <- x$co
        else if ("c1" %in% att)
            X <- x$c1
        else stop("Can't find scores")
    }
    else if (is.numeric(x)) {
        X <- as.matrix(x)
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
    if (is.null(rownames(X))) {
        root <- substr(display, 1, 4)
        rownames(X) <- paste(root, 1:nrow(X), sep = "")
    }
    if (is.null(colnames(X))) 
        colnames(X) <- paste("Dim", 1:ncol(X), sep = "")
    if (!missing(choices)) {
        choices <- choices[choices <= ncol(X)]
        X <- X[, choices, drop = FALSE]
    }
    X <- as.matrix(X)
    X
}
