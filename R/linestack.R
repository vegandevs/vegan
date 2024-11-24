`linestack` <-
    function (x, labels, cex = 0.8, side = "right", hoff = 2, air = 1.1,
              at = 0, add = FALSE, axis = FALSE, ...)
{
    if (length(at) > 1 || length(hoff) > 1 || length(air) > 1 || length(cex) > 1)
        stop("only one value accepted for arguments 'cex', 'hoff', 'air' and 'at'")
    x <- drop(x)
    n <- length(x)
    misslab <- missing(labels)
    if (misslab) {
        labels <- names(x)
    }
    if (!is.expression(labels) && !is.character(labels)) {
        labels <- as.character(labels)  # coerce to character only if not expressions
    }
    nlab <- length(labels)
    if (!misslab && n != nlab) {
        stop(gettextf(
            "wrong number of supplied 'labels: expected %d, got %d", n, nlab))
    }
    side <- match.arg(side, c("right", "left"))
    op <- par(xpd = TRUE)
    on.exit(par(op))
    ord <- order(x)
    x <- x[ord]
    labels <- labels[ord]
    pos <- numeric(n)
    if (!add) {
        plot(pos, x, type = "n", axes = FALSE, xlab = "", ylab = "", ...)
    }
    hoff <- hoff * strwidth("m")
    ht <- air * strheight(labels, cex = cex)
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
    ## plot text in the original order for correct matching of vectors
    ## of graphical parameters (...).
    o <- order(ord)
    if (side == "right") {
        text(at + hoff, pos[o], labels[o], pos = 4, cex = cex, offset = 0.2,
             ...)
        segments(at, x, at + hoff, pos)
    }
    else if (side == "left") {
        text(at - hoff, pos[o], labels[o], pos = 2, cex = cex, offset = 0.2,
             ...)
        segments(at, x, at - hoff, pos)
    }
    if (axis)
        axis(if (side == "right")
             2
        else 4, pos = at, las = 2)
    invisible(pos[o])
}
