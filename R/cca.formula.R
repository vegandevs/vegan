"cca.formula" <-
function (formula, data, ...) 
{
    if (missing(data)) {
        data <- parent.frame()
    } else {
        data <- ordiGetData(match.call(), environment(formula))
    }
    d <- ordiParseFormula(formula, data)
    sol <- cca.default(d$X, d$Y, d$Z)
    if (!is.null(sol$CCA)) 
        sol$CCA$centroids <- centroids.cca(sol$CCA$wa, d$modelframe, 
            sol$rowsum)
    if (!is.null(sol$CCA$alias)) 
        sol$CCA$centroids <- unique(sol$CCA$centroids)
    if (!is.null(sol$CCA$centroids)) {
        rs <- rowSums(sol$CCA$centroids^2)
        sol$CCA$centroids <- sol$CCA$centroids[rs > 1e-04, , 
            drop = FALSE]
    }
    sol$terms <- d$terms
    sol$terminfo <- ordiTerminfo(d, data)
    sol$call <- match.call()
    sol$call[[1]] <- as.name("cca")
    sol$call$formula <- formula(d$terms, width.cutoff = 500)
    sol
}
