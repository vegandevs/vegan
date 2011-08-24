`predict.specaccum` <-
    function(object, newdata, interpolation = c("linear", "spline"), ...)
{
    if (missing(newdata))
        out <- object$richness
    else {
        interpolation <- match.arg(interpolation)
        newdata <- drop(as.matrix(newdata))
        if (length(dim(newdata)) > 1)
            stop("function accepts only one variable as 'newdata'")
        if (interpolation == "linear")
            out <- approx(x = object$sites, y = object$richness,
                          xout = newdata, rule = 1)$y
        else
            out <- spline(x = object$sites, y = object$richness,
                          xout = newdata, ...)$y
    }
    out
}
