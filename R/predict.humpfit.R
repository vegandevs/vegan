"predict.humpfit" <-
    function(object, newdata = NULL, ...)
{
    if (is.null(newdata))
        return(fitted(object))
    else {
        p <- coef(object)
        x <- unlist(newdata)
        x <- ifelse(x < p[1], x/p[1], p[1]*p[1]/x/x)
        fv <- p[3]*log(1 + p[2]*x/p[3])
    }
    fv
}
