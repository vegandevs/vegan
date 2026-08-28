`text.decorana` <-
    function (x, display = c("sites", "species"), labels, choices = 1:2,
              origin = TRUE, select, ...)
{
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
    text.ordiplot(x, what = display, ...)
}
