`points.cca` <-
    function (x, display = "sites", choices = c(1, 2), scaling = "species",
              arrow.mul, head.arrow = 0.05, select, const,
              correlation = FALSE, hill = FALSE, ...)
{
    if (length(display) > 1)
        stop("only one 'display' item can be added in one command")
    pts <- scores(x, choices = choices, display = display, scaling = scaling,
                  const, correlation = correlation, hill = hill, tidy=FALSE,
                  droplist = FALSE)
    class(pts) <- "ordiplot"
    if (!missing(select))
        pts[[1]] <- .checkSelect(select, pts[[1]])
    points.ordiplot(pts, what = names(pts), ...)
    invisible()
}

