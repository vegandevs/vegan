"fitted.rda" <-
    function (object, model = c("CCA", "CA"), type = c("response", "working"), ...) 
{
    type <- match.arg(type)
    model <- match.arg(model)
    Xbar <- object[[model]]$Xbar
    if (model == "CCA") 
        Xbar <- qr.fitted(object$CCA$QR, Xbar)
    if (type == "response") {
        cent <- attr(Xbar, "scaled:center")
        scal <- attr(Xbar, "scaled:scale")
        if (!is.null(scal)) {
            Xbar <- sweep(Xbar, 2, scal, "*")
            attr(Xbar, "scaled:scale") <- NULL
        }
        Xbar <- sweep(Xbar, 2, cent, "+")
        attr(Xbar, "scaled:center") <- NULL
    } else {
        Xbar <- Xbar/sqrt(nrow(Xbar)-1)
    }
    Xbar
}
