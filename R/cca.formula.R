`cca.formula` <-
    function (formula, data, na.action = na.fail, subset = NULL, ...)
{
    if (missing(data)) {
        data <- parent.frame()
    } else {
        data <- eval(match.call()$data, environment(formula),
                     enclos = .GlobalEnv)
    }
    d <- ordiParseFormula(formula, data = data, na.action = na.action,
                          subset = substitute(subset))
    sol <- cca.default(d$X, d$Y, d$Z)
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
        centroids <- centroids.cca(sol$CCA$wa, d$modelframe,
                                   sol$rowsum)
        if (!is.null(sol$CCA$alias))
            centroids <- unique(centroids)
        ## See that there really are centroids
        if (!is.null(centroids)) {
            rs <- rowSums(centroids^2)
            centroids <- centroids[rs > 1e-04,, drop = FALSE]
            if (length(centroids) == 0)
                centroids <- NULL
        }
        if (!is.null(centroids))
            sol$CCA$centroids <- centroids
    }
    ## replace cca.default call
    call <- match.call()
    call[[1]] <- as.name("cca")
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
    class(sol) <- "cca"
    sol
}
