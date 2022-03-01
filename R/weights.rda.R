`weights.rda` <-
    function (object, display = "sites", ...)
{
    display <- match.arg(display, c("sites", "species", "lc",
                                    "wa"))
    if (display %in% c("sites", "lc", "wa")) {
        n <- nobs(object)
        if (!is.null(object$na.action) &&
            inherits(object$na.action, "exclude"))
            n <- n + length(object$na.action)
    }
    else n <- length(object$colsum)
    rep(1, n)
}
