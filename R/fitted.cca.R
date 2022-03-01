`fitted.cca` <-
    function (object, model = c("CCA","CA","pCCA"),
              type = c("response", "working"), ...)
{
    type <- match.arg(type)
    model <- match.arg(model)
    if (is.null(object[[model]]))
        stop(gettextf("component '%s' does not exist", model))
    Xbar <- ordiYbar(object, model)
    if (type == "response") {
        gtot <- object$grand.total
        rc <- object$rowsum %o% object$colsum
        Xbar <- (Xbar * sqrt(rc) + rc) * gtot
    }
    Xbar
}
