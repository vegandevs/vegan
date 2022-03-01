`fitted.rda` <-
    function (object, model = c("CCA", "CA", "pCCA"),
              type = c("response", "working"), ...)
{
    type <- match.arg(type)
    model <- match.arg(model)
    if (is.null(object[[model]]))
        stop(gettextf("component '%s' does not exist", model))
    Xbar <- ordiYbar(object, model)
    if (type == "response") {
        cent <- attr(Xbar, "scaled:center")
        scal <- attr(Xbar, "scaled:scale")
        if (!is.null(scal)) {
            Xbar <- sweep(Xbar, 2, scal, "*")
            attr(Xbar, "scaled:scale") <- NULL
        }
        Xbar <- Xbar * sqrt(nrow(Xbar) - 1)
        Xbar <- sweep(Xbar, 2, cent, "+")
        attr(Xbar, "scaled:center") <- NULL
    }
    Xbar
}
