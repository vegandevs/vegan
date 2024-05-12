"summary.decorana" <-
    function (object, digits = 3, origin = TRUE, display = c("both",
                                                 "species", "sites", "none"), ...)
{
    warning("'summary' is deprecated: use 'scores' for scores, 'weights' for weights")
    display <- match.arg(display)
    print(object)
    if (origin) {
        object$cproj <- sweep(object$cproj, 2, object$origin,
                              "-")
        object$rproj <- sweep(object$rproj, 2, object$origin,
                              "-")
    }
    tmp <- list()
    if (display == "both" || display == "species") {
        tmp$spec.scores <- object$cproj
        tmp$spec.priorweights <- object$v
        tmp$spec.totals <- object$adotj
    }
    if (display == "both" || display == "sites") {
        tmp$site.scores <- object$rproj
        tmp$site.totals <- object$aidot
    }
    tmp$digits <- digits
    class(tmp) <- "summary.decorana"
    tmp
}
