`hiersimu.formula` <-
    function(formula, data, FUN, location = c("mean", "median"),
             relative = FALSE, drop.highest = FALSE, nsimul=99, ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    lhs <- tmp$lhs
    rhs <- tmp$rhs

    ## run simulations
    sim <- hiersimu.default(lhs, rhs, FUN = FUN, location = location,
                            relative = relative, drop.highest = drop.highest,
                            nsimul = nsimul, ...)
    call <- match.call()
    call[[1]] <- as.name("hiersimu")
    attr(sim, "call") <- call
    sim
}

