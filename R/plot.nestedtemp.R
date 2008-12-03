`plot.nestedtemp` <-
    function (x, kind = c("temperature", "incidendce"),
              col = rev(heat.colors(100)), names = FALSE,  
              ...) 
{
    kind <- match.arg(kind)
    if (kind == "temperature") 
        z <- x$u
    else z <- x$comm
    z <- t(z[nrow(z):1, ])
    image(z, axes = FALSE, col = col, ...)
    box()
    lines(x$smooth$x, 1 - x$smooth$y)
    if (length(names) == 1)
        names <- rep(names, 2)
    if (names[1]) {
        axis(2, at = seq(1, 0, len = nrow(x$u)), labels = rownames(x$u),
            las = 2, ...)
    }
    if (names[2]) {
        axis(3, at = seq(0, 1, len = ncol(x$u)), labels = colnames(x$u),
            las = 2, ...)
    }
}
