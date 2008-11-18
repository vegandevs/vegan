"text.ordiplot" <-
    function (x, what, labels, select, ...) 
{
    x <- scores(x, what)
    if (!missing(labels))
        rownames(x) <- labels
    if (!missing(select)) 
        x <- x[select, , drop = FALSE]
    text(x, labels = rownames(x), ...)
    invisible()
}
