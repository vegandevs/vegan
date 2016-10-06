`envfit.formula` <-
    function(formula, data, ...)
{
    if (missing(data))
        data <- environment(formula)
    X <- formula[[2]]
    X <- eval(X, environment(formula), enclos = .GlobalEnv)
    formula[[2]] <- NULL
    P <- model.frame(formula, data, na.action = na.pass)
    envfit.default(X, P, ...)
}
