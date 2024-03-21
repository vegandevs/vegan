`scores.orditkplot` <-
    function(x, display, ...)
{
    .Deprecated(msg = "function was moved to CRAN package vegan3d")
    if (!missing(display) && !is.na(pmatch(display, "labels")))
        x$labels
    else
        x$points
}

