`bioenv.formula` <-
    function (formula, data, ...) 
{
    if (missing(data)) 
        data <- parent.frame()
    fla <- formula
    comm <- formula[[2]]
    comm <- eval(comm, data, parent.frame())
    formula[[2]] <- NULL
    mf <- model.frame(formula, data, na.action = NULL)
    if (any(sapply(mf, function(x) is.factor(x) || !is.numeric(x)))) 
        stop("bioenv applies only to numeric variables")
    env <- attr(mf, "terms")
    attr(env, "intercept") <- 0
    env <- model.matrix(env, mf)
    out <- bioenv(comm, env, ...)
    out$formula <- fla
    out$call <- match.call()
    out$call[[1]] <- as.name("bioenv")
    out
}
