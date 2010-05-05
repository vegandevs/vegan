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
    if (is.data.frame(env)) {
        facts <- unlist(lapply(env, is.factor))
        if (sum(facts)) {
            Pfac <- env[, facts, drop = FALSE]
            P <- env[, !facts, drop = FALSE]
            if (length(P)) {
                if (permutations) {
                    if (!exists(".Random.seed", envir = .GlobalEnv, 
                                inherits = FALSE)) {
                        runif(1)
                    }
                    seed <- get(".Random.seed", envir = .GlobalEnv, 
                                inherits = FALSE)
                }
                vectors <- vectorfit(X, P, permutations, strata, 
                                     choices, w = w, ...)
            }
            if (!is.null(seed)) {
                assign(".Random.seed", seed, envir = .GlobalEnv)
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
