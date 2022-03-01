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
    ## take care that the first column is a unique identifier for rows
    ## and the last column is constant for pooling all rows together
    if (length(unique(rhs[,1])) < nrow(rhs))
        rhs <- cbind("unit" = factor(seq_len(nrow(rhs))), rhs)
    if (length(unique(rhs[, ncol(rhs)])) > 1)
        rhs <- cbind(rhs, "all" = factor(rep(1, nrow(rhs))))
    attr(rhs, "terms") <- NULL
    list(lhs=lhs, rhs=rhs)
}

