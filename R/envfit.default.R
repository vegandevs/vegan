"envfit.default" <-
    function (ord, env, permutations = 999, strata, choices = c(1, 2), 
             display = "sites", w = weights(ord), na.rm = FALSE, ...) 
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
        env <- env[keep,, drop=FALSE]
        na.action <- structure(seq_along(keep)[!keep], class="omit")
    }
    ## make permutation matrix for all variables handled in the next loop
    nr <- nrow(X)
    if (length(permutations) == 1) {
        if (permutations > 0 ) {
            arg <- if (missing(strata)) NULL else strata
            permutations <- t(replicate(permutations,
                                  permuted.index(nr, strata=arg)))
        }
    } else {
        permat <- as.matrix(permutations)
        if (ncol(permat) != nr)
            stop(gettextf("'permutations' have %d columns, but data have %d rows",
                          ncol(permat), nr))
    }
    if (is.data.frame(env)) {
        facts <- sapply(env, is.factor)
        if (sum(facts)) {  # have factors
            Pfac <- env[, facts, drop = FALSE]
            P <- env[, !facts, drop = FALSE]
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
