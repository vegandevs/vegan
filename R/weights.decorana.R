"weights.decorana" <-
    function(object, display="sites", ...)
{
    display <- match.arg(display, c("sites","species"))
    if (display == "sites") object$aidot
    else object$adotj
}
