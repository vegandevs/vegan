`adipart.formula` <-
    function(formula, data, index=c("richness", "shannon", "simpson"),
             weights=c("unif", "prop"), relative = FALSE, nsimul=99,
             method = "r2dtable", ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    ## run simulations
    sim <- adipart.default(tmp$lhs, tmp$rhs, index = index, weights = weights,
                           relative = relative, nsimul = nsimul,
                           method = method, ...)
    call <- match.call()
    call[[1]] <- as.name("adipart")
    attr(sim, "call") <- call
    sim
}
