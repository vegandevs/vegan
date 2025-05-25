`text.metaMDS` <-
    function (x, display = c("sites", "species"), labels,
              choices = c(1, 2), shrink = FALSE, select, cex = 0.7, ...)
{
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, shrink = shrink)
    if (!missing(select))
        x <- .checkSelect(select, x)
    if (!missing(labels))
        rownames(x) <- labels
    text.ordiplot(x, what = display, labels = rownames(x), cex = cex, ...)
    invisible()
}
