`text.orditkplot` <-
    function(x, ...)
{
    text(x$labels, label = rownames(x$labels), ...)
}

