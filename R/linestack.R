"linestack" <-
    function (x, labels, cex = 0.8, side = "right", hoff = 2, air = 1.1, 
              at = 0, add = FALSE, axis = FALSE, ...) 
{
    if (!missing(labels) && length(labels == 1) && pmatch(labels, 
                                   c("right", "left"), nomatch = FALSE)) {
        side <- labels
        labels <- NULL
        warning("argument 'label' is deprecated: use 'side'")
    }
    side <- match.arg(side, c("right", "left"))
    x <- drop(x)
    if (!missing(labels) && !is.null(labels)) 
        names(x) <- labels
    else if (is.null(names(x)))
        names(x) <- rep("", length(x))
    op <- par(xpd = TRUE)
    ord <- order(x)
    x <- x[ord]
    n <- length(x)
    pos <- numeric(n)
    if (!add) {
        plot(pos, x, type = "n", axes = FALSE, xlab = "", ylab = "", ...)
    }
    hoff <- hoff * strwidth("m")
    ht <- air * strheight(names(x), cex = cex)
    mid <- (n + 1)%/%2
    pos[mid] <- x[mid]
    if (n > 1) {
        for (i in (mid + 1):n) {
            pos[i] <- max(x[i], pos[i - 1] + ht[i])
        }
    }
    if (n > 2) {
        for (i in (mid - 1):1) {
            pos[i] <- min(x[i], pos[i + 1] - ht[i])
        }
    }
    segments(at, x[1], at, x[n])
    if (side == "right") {
        text(at + hoff, pos, names(x), pos = 4, cex = cex, offset = 0.2, 
             ...)
        segments(at, x, at + hoff, pos)
    }
    else if (side == "left") {
        text(at - hoff, pos, names(x), pos = 2, cex = cex, offset = 0.2, 
             ...)
        segments(at, x, at - hoff, pos)
    }
    if (axis) 
        axis(if (side == "right") 
             2
        else 4, pos = at, las = 2)
    par(op)
    invisible(pos[order(ord)])
}
