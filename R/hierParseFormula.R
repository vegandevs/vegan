`hierParseFormula` <-
    function (formula, data)
{
    lhs <- formula[[2]]
    if (any(attr(terms(formula, data = data), "order") > 1))
        stop("interactions are not allowed")
    lhs <- as.matrix(eval(lhs, environment(formula), parent.frame()))
    formula[[2]] <- NULL
    rhs <- model.frame(formula, data, drop.unused.levels = TRUE)
    rhs[] <- lapply(rhs, function(u) {
        if (!is.factor(u))
            u <- factor(u)
        u
    })
    if (length(rhs) < 2)
        stop("at least 2 hierarchy levels are needed")
    attr(rhs, "terms") <- NULL
    list(lhs=lhs, rhs=rhs)
}

