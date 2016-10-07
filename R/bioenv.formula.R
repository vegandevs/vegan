`bioenv.formula` <-
    function (formula, data, ...)
{
    if (missing(data))
        data <- environment(formula)
    fla <- formula
    comm <- formula[[2]]
    comm <- eval(comm, environment(formula), parent.frame())
    formula[[2]] <- NULL
    env <- model.frame(formula, data, na.action = NULL)
    out <- bioenv(comm, env, ...)
    out$formula <- fla
    out$call <- match.call()
    out$call[[1]] <- as.name("bioenv")
    out
}
