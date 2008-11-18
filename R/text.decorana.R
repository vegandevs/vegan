"text.decorana" <-
    function (x, display = c("sites", "species"), labels, choices = 1:2,
              origin = TRUE, select, ...)
{
    localText <- function(..., shrink, origin, scaling, triangular)
        segments(...)
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, origin = origin,
                ...)
    if (!missing(labels))
        rownames(x) <- labels
    if (!missing(select))
        x <- x[select, , drop = FALSE]
    localText(x, rownames(x), ...)
    invisible()
}
