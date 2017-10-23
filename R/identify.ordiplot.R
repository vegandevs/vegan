`identify.ordiplot` <-
    function (x, what, labels, ...)
{
    x <- scores(x, what)
    if (missing(labels))
        labels <- rownames(x)
    identify(x, labels = labels, ...)
}
