"envfit.default" <-
    function (X, P, permutations = 0, strata, choices = c(1, 2), 
              ...) 
{
    vectors <- NULL
    factors <- NULL
    seed <- NULL
    if (is.data.frame(P)) {
        facts <- unlist(lapply(P, is.factor))
        if (sum(facts)) {
            Pfac <- P[, facts, drop = FALSE]
            P <- P[, !facts, drop = FALSE]
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
                                     choices, ...)
            }
            if (!is.null(seed)) {
                assign(".Random.seed", seed, envir = .GlobalEnv)
            }
            factors <- factorfit(X, Pfac, permutations, strata, 
                                 choices, ...)
            sol <- list(vector = vectors, factors = factors)
        }
        else vectors <- vectorfit(X, P, permutations, strata, 
                                  choices, ...)
    }
    else vectors <- vectorfit(X, P, permutations, strata, choices)
    sol <- list(vectors = vectors, factors = factors)
    class(sol) <- "envfit"
    sol
}
