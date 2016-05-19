`hiersimu.formula` <-
    function(formula, data, FUN, location = c("mean", "median"),
             relative = FALSE, drop.highest = FALSE, nsimul=99,
             method = "r2dtable", ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    ## run simulations
    sim <- hiersimu.default(tmp$lhs, tmp$rhs, FUN = FUN, location = location,
                            relative = relative, drop.highest = drop.highest,
                            nsimul = nsimul, method = method, ...)
    call <- match.call()
    call[[1]] <- as.name("hiersimu")
    attr(sim, "call") <- call
    sim
}

