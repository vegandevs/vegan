"text.decorana" <-
    function (x, display = c("sites", "species"), labels, choices = 1:2,
              origin = TRUE, select, ...)
{
    localText <- function(..., shrink, origin, scaling, triangular)
        text(...)
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, origin = origin,
                ...)
    if (!missing(labels))
        rownames(x) <- labels
    if (!missing(select))
        x <- .checkSelect(select, x)
    localText(x, rownames(x), ...)
    invisible()
}
