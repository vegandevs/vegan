`text.orditkplot` <-
    function(x, ...)
{
    text(x$labels, labels = rownames(x$labels), ...)
}

