`ordilabel` <-
    function(x, display, labels, choices = c(1,2), priority,
             cex = 0.8, fill = "white", border = NULL,  ...)
{
    x <- scores(x, display = display, choices = choices, ...)
    if (missing(labels))
        labels <- rownames(x)
    if (!missing(priority)) {
        ord <- order(priority)
        x <- x[ord, ]
        labels <- labels[ord]
    }
    em <- strwidth("m", cex = cex, ...)
    ex <- strheight("x", cex = cex, ...)
    w <- (strwidth(labels, cex=cex,...) + em/1.5)/2
    h <- (strheight(labels, cex = cex, ...) + ex/1.5)/2
    for (i in 1:nrow(x)) {
        polygon(x[i,1] + c(-1,1,1,-1)*w[i], x[i,2] + c(-1,-1,1,1)*h[i],
                col = fill, border = border)
        text(x[i,1], x[i,2], labels = labels[i], cex = cex, ...)
    }
    invisible(x)
}

