"lines.radline" <-
    function (x, ...) 
{
    lin <- fitted(x)
    rnk <- seq(along = lin)
    lines(rnk, lin, ...)
    invisible()
}
