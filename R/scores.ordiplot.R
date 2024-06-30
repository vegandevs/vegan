`scores.ordiplot`  <-
    function (x, display = "sites", ...)
{
    if (length(x) == 1) {
        attr(x[[1]], "score") <- names(x)
        return(x[[1]])
    }
    items <- names(x)
    items <- items[!is.na(items)]
    display <- match.arg(display, items)
    x <- x[[display]]
    attr(x, "score") <- display
    x
}
