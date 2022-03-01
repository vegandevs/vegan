`multipart.formula` <-
    function(formula, data, index=c("renyi", "tsallis"), scales = 1,
             global = FALSE, relative = FALSE, nsimul=99,
             method = "r2dtable", ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    ## run simulations
    sim <- multipart.default(tmp$lhs, tmp$rhs, index = index, scales = scales,
                             global = global, relative = relative,
                             nsimul = nsimul, method = method, ...)
    call <- match.call()
    call[[1]] <- as.name("multipart")
    attr(sim, "call") <- call
    sim
}
