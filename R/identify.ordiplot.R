`identify.ordiplot` <-
    function (x, what, labels, ...)
{
    x <- scores(x, display = what)
    if (missing(labels))
        labels <- rownames(x)
    identify(x, labels = labels, ...)
}
