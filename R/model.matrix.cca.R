`model.matrix.cca` <-
    function (object, ...) 
{
    if (inherits(object, "prc"))
        stop("model.matrix does not work with 'prc' results")
    call <- object$call
    m <- match(c("formula", "data", "na.action", "subset"), names(call), 
        0)
    call <- call[c(1, m)]
    call[[1]] <- as.name("ordiParseFormula")
    out <- eval(call, environment(), parent.frame())[c("Z", "Y")]
    m <- list()
    if (!is.null(out$Z))
        m$Conditions <- out$Z
    if (!is.null(out$Y))
        m$Constraints <- out$Y
    if (length(m) == 1)
        m <- m[[1]]
    m
}
