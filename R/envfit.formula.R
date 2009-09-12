"envfit.formula" <-
    function(formula, data, ...)
{
    if (missing(data))
        data <- parent.frame()
    X <- formula[[2]]
    X <- eval(X, data, parent.frame())
    formula[[2]] <- NULL
    P <- model.frame(formula, data, na.action = na.pass)
    envfit(X, P, ...)
}
