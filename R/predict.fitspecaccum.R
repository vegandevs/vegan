`predict.fitspecaccum` <-
    function(object, newdata, ...)
{
    mods <- object$models
    if (!missing(newdata)) {
        newdata <- drop(as.matrix(newdata))
        if (length(dim(newdata)) > 1)
            stop("function accepts only one variable as 'newdata'")
        drop(sapply(mods, predict, newdata = data.frame(x = newdata), ...))
    } else {
        drop(sapply(mods, predict, ...))
    }
}
