`lines.radline` <-
    function (x, ...) 
{
    lin <- fitted(x)
    rnk <- seq(along = lin)
    lines(rnk, lin, ...)
    invisible()
}

`lines.radfit` <-
    function(x, ...)
{
    matlines(fitted(x), ...)
}
