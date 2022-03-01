"text.metaMDS" <-
    function (x, display = c("sites", "species"), labels,
              choices = c(1, 2), shrink = FALSE, select, ...)
{
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, shrink = shrink)
    if (!missing(labels))
        rownames(x) <- labels
    if (!missing(select))
        x <- .checkSelect(select, x)
    text(x, labels = rownames(x), ...)
    invisible()
}
