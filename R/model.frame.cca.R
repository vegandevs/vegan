`model.frame.cca` <-
    function (formula, ...) 
{
    call <- formula$call
    m <- match(c("formula", "data", "na.action", "subset"), names(call), 
        0)
    call <- call[c(1, m)]
    call[[1]] <- as.name("ordiParseFormula")
    out <- eval(call, parent.frame())
    mf <- out$modelframe
    attr(mf, "terms") <- out$terms.expand
    if (!is.null(out$na.action)) 
        attr(mf, "na.action") <- out$na.action
    mf
}
