`text.decorana` <-
    function (x, display = c("sites", "species"), labels, choices = 1:2,
              origin = TRUE, select, ...)
{
    ## do we need localText???
    localText <- function(..., shrink, origin, scaling, triangular)
        text.ordiplot(...)
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, origin = origin,
                ...)
    ## text.ordiplot expects a list of one matrix
    x <- list(x)
    names(x) <- display
    if (!missing(select))
        x <- .checkSelect(select, x)
    if (!missing(labels))
        rownames(x) <- labels
    localText(x, what = display, ...)
    invisible()
}
