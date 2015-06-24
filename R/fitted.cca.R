`fitted.cca` <-
    function (object, model = c("CCA","CA","pCCA"), type = c("response", "working"),
              ...) 
{
    type <- match.arg(type)
    model <- match.arg(model)
    if (is.null(object[[model]]))
        stop("component ", model, " does not exist")
    gtot <- object$grand.total
    rc <- object$rowsum %o% object$colsum
    if (model == "pCCA")
        Xbar <- object$pCCA$Fit
    else
        Xbar <- object[[model]]$Xbar
    if (model == "CCA")
        Xbar <- qr.fitted(object$CCA$QR, Xbar)
    if (type == "response")
        Xbar <- (Xbar * sqrt(rc)   + rc) * gtot 
    Xbar
}
