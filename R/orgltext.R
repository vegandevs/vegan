"orgltext" <-
    function (object, text, display = "sites", choices = 1:3, justify = "center",  
              adj = 0.5, ...) 
{
    x <- scores(object, display = display, choices = choices, 
                ...)
    if (missing(text)) 
        text <- rownames(x)
    rgl.texts(x[, 1], x[, 2], x[, 3], text, adj = adj,  ...)
    invisible()
}
