`envfit.default` <-
    function (ord, env, permutations = 999, strata = NULL, choices = c(1, 2),
             display = "sites", w = weights(ord, display), na.rm = FALSE, ...)
{
    weights.default <- function(object, ...) NULL
    w < eval(w)
    vectors <- NULL
    factors <- NULL
    seed <- NULL
    X <- scores(ord, display = display, choices = choices, ...)
    keep <- complete.cases(X) & complete.cases(env)
    if (any(!keep)) {
        if (!na.rm)
            stop("missing values in data: consider na.rm = TRUE")
        X <- X[keep,, drop=FALSE]
        ## drop any lost levels, explicitly don't include NA as a level
        env <- droplevels(env[keep,, drop=FALSE], exclude = NA)
        na.action <- structure(seq_along(keep)[!keep], class="omit")
    }
    ## make permutation matrix for all variables handled in the next loop
    nr <- nrow(X)
    permat <-  getPermuteMatrix(permutations, nr, strata = strata)
    if (ncol(permat) != nr)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                      ncol(permat), nr))

    if (is.data.frame(env)) {
        vects <- sapply(env, is.numeric)
        if (any(!vects)) {  # have factors
            Pfac <- env[, !vects, drop = FALSE]
            P <- env[, vects, drop = FALSE]
            if (length(P)) { # also have vectors
                vectors <- vectorfit(X, P, permutations, strata,
                                     choices, w = w, ...)
            }
            factors <- factorfit(X, Pfac, permutations, strata,
                                         choices, w = w, ...)
            sol <- list(vector = vectors, factors = factors)
            }
        else vectors <- vectorfit(X, env, permutations, strata,
                                  choices, w = w, ...)
    }
    else vectors <- vectorfit(X, env, permutations, strata,
                              choices, w = w, ...)
    sol <- list(vectors = vectors, factors = factors)
    if (!is.null(na.action))
        sol$na.action <- na.action
    class(sol) <- "envfit"
    sol
}
