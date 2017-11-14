`text.ordiplot`  <-
    function (x, what, labels, select, ...)
{
    sco <- scores(x, what)
    if (!missing(labels))
        rownames(sco) <- labels
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    text(sco, labels = rownames(sco), ...)
    invisible(x)
}
