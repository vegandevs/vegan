`adipart.formula` <-
    function(formula, data, index=c("richness", "shannon", "simpson"),
             weights=c("unif", "prop"), relative = FALSE, nsimul=99, ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    lhs <- tmp$lhs
    rhs <- tmp$rhs

    ## run simulations
    sim <- adipart.default(lhs, rhs, index = index, weights = weights,
                           relative = relative, nsimul = nsimul, ...)
    call <- match.call()
    call[[1]] <- as.name("adipart")
    attr(sim, "call") <- call
    sim
}
