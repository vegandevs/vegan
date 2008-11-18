"orglpoints" <-
    function (object, display = "sites", choices = 1:3, ...) 
{
    x <- scores(object, display = display, choices = choices, ...)
    rgl.points(x[,1], x[,2], x[,3], ...)
    invisible()
}

