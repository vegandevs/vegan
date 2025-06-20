`wrda` <-
    function (formula, data, w, na.action = na.fail,
              subset = NULL, ...)
{
    if (missing(data)) {
        data <- parent.frame()
    } else {
        data <- eval(match.call()$data, environment(formula),
                     enclos = .GlobalEnv)
    }
    ## Use equal weights with the same result as in rda if w not given
    if (missing(w))
        w <- rep(1/(nrow(data)-1), nrow(data))
    d <- ordiParseFormula(formula, data = data, na.action = na.action,
                          subset = substitute(subset))
    sol <- ordConstrained(d$X, d$Y, d$Z, "wrda", w)
    sol$inertia <- "weighted variance"
    sol$CCA$centroids <- getCentroids(sol, d$modelframe)
    ## replace rda.default call
    call <- match.call()
    ## call[[1]] <- as.name("rda")
    call$formula <- formula(d$terms)
    sol$call <- call
    if (!is.null(d$na.action)) {
        sol$na.action <- d$na.action
        sol <- ordiNAexclude(sol, d$excluded)
    }
    if (!is.null(d$subset))
        sol$subset <- d$subset
    ## drops class in c()
    sol <- c(sol,
             list(terms = d$terms,
                  terminfo = ordiTerminfo(d, d$modelframe)))
    class(sol) <- c("wrda", "rda", "cca")
    sol
}
