"weights.cca" <-
    function (object, display = "sites", ...) 
{
    display <- match.arg(display, c("sites", "species", "lc", "wa"))
    if (display %in% c("sites", "lc", "wa")) 
        object$rowsum
    else object$colsum
}
