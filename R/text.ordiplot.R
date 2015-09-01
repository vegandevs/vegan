"text.ordiplot" <-
    function (x, what, labels, select, ...)
{
    x <- scores(x, what)
    if (!missing(labels))
        rownames(x) <- labels
    if (!missing(select))
        x <- .checkSelect(select, x)
    text(x, labels = rownames(x), ...)
    invisible()
}
