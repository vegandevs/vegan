`scores.orditkplot` <-
    function(x, display, ...)
{
    if (!missing(display) && !is.na(pmatch(display, "labels")))
        x$labels
    else
        x$points
}

