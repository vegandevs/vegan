"scores.ordiplot" <-
function (x, display = "sites", ...) 
{
    if (length(x) == 1)
        return(x[[1]])
    items <- names(x)
    items <- items[!is.na(items)]
    display <- match.arg(display, items)
    cmd <- paste("x", display, sep = "$")
    eval(parse(text = cmd))
}
