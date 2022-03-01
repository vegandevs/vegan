"weights.cca" <-
    function (object, display = "sites", ...) 
{
    display <- match.arg(display, c("sites", "species", "lc", "wa"))
    if (display %in% c("sites", "lc", "wa")) {
        if (!is.null(object$na.action) &&
            inherits(object$na.action, "exclude")) {
            object$rowsum <- napredict(object$na.action, object$rowsum)
            object$rowsum[object$na.action] <- object$rowsum.excluded
            }
        object$rowsum
    }
    else object$colsum
}
