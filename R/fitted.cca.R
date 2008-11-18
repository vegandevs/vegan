"fitted.cca" <-
function (object, model = c("CCA","CA"), type = c("response", "working"), ...) 
{
    type <- match.arg(type)
    model <- match.arg(model)
    gtot <- object$grand.total
    rc <- object$rowsum %o% object$colsum
    Xbar <- object[[model]]$Xbar
    if (model == "CCA")
        Xbar <- qr.fitted(object$CCA$QR, Xbar)
    if (type == "response")
        Xbar <- (Xbar * sqrt(rc)   + rc) * gtot 
    Xbar
}
