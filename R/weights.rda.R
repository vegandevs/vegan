"weights.rda" <-
    function (object, display = "sites", ...) 
{
    display <- match.arg(display, c("sites", "species", "lc", 
                                    "wa"))
    n <- if (display %in% c("sites", "lc", "wa")) 
        max(nrow(object$CA$Xbar), nrow(object$CCA$Xbar))
    else max(ncol(object$CA$Xbar), ncol(object$CCA$Xbar))
    rep(1, n)
}
