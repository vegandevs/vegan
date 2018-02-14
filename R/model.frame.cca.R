`model.frame.cca` <-
    function (formula, ...)
{
    call <- formula$call
    m <- match(c("formula", "data", "na.action", "subset"), names(call),
        0)
    call <- call[c(1, m)]
    ## did we succeed? Fails if we have no formula, in prc and if
    ## there was no data= argument
    if (is.null(call$data))
        stop("no sufficient information to reconstruct model frame")
    ## subset must be evaluated before ordiParseFormula
    if (!is.null(call$subset))
        call$subset <- formula$subset
    if (is.null(call$na.action))
        call$na.action <- na.pass
    data <- eval(call$data, environment(call$formula), .GlobalEnv)
    out <- ordiParseFormula(call$formula, data, na.action = call$na.action,
                            subset = call$subset)
    mf <- out$modelframe
    attr(mf, "terms") <- out$terms.expand
    if (!is.null(out$na.action))
        attr(mf, "na.action") <- out$na.action
    mf
}
