"weights.rda" <-
    function (object, display = "sites", ...) 
{
    display <- match.arg(display, c("sites", "species", "lc", 
                                    "wa"))
    if (display %in% c("sites", "lc", "wa")) {
        n <- max(nrow(object$CA$Xbar), nrow(object$CCA$Xbar))
        if (!is.null(object$na.action) &&
            inherits(object$na.action, "exclude"))
            n <- n + length(object$na.action)
    }
    else n <- max(ncol(object$CA$Xbar), ncol(object$CCA$Xbar))
    rep(1, n)
}
