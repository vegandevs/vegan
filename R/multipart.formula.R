`multipart.formula` <-
    function(formula, data, index=c("renyi", "tsallis"), scales = 1,
             global = FALSE, relative = FALSE, nsimul=99, ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    lhs <- tmp$lhs
    rhs <- tmp$rhs

    ## run simulations
    sim <- multipart.default(lhs, rhs, index = index, scales = scales,
                             global = global, relative = relative,
                             nsimul = nsimul, ...)
    call <- match.call()
    call[[1]] <- as.name("multipart")
    attr(sim, "call") <- call
    sim
}
