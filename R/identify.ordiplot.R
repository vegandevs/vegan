"identify.ordiplot" <-
    function (x, what, labels, ...) 
{
    x <- scores(x, what)
    if (missing(labels))
        labels <- rownames(x)
    out <- identify(x, labels = labels, ...)
    out
}
