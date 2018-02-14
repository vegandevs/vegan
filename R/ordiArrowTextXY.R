### Location of the text at the point of the arrow. 'x' are the
### coordinates of the arrow heads, and 'labels' are the text used to
### label these heads, '...' passes arguments (such as 'cex') to
### strwidth() and strheight().
`ordiArrowTextXY` <- function (x, labels, display, choices = c(1,2),
                               rescale = TRUE, fill = 0.75, ...) {
    ## handle x, which we try with scores, but also retain past usage of
    ## a two column matrix
    X <- if (is.matrix(x)) {
        nc <- NCOL(x)
        if (nc != 2L) {
            stop("a two-column matrix of coordinates is required")
        }
        x
    } else {
        if (inherits(x, "envfit")) {
            scores(x, display = "vectors", ...)[, 1:2]
        } else {
            scores(x, display = display, choices = choices, ...)
        }
        if (!rescale) {
            warning("extracted scores usually need rescaling but you set 'rescale = FALSE' - \nconsider using 'rescale = TRUE', the default")
        }
    }

    ## find multiplier to fill if rescaling
    if (rescale) {
        mul <- ordiArrowMul(X, fill = fill)
        X <- X * mul
    }

    if (missing(labels)) {
        rnames <- rownames(X)
        labels <- if (is.null(rnames)) {
            paste("V", seq_len(NROW(X)))
        } else {
            rnames
        }
    }

    w <- strwidth(labels, ...)
    h <- strheight(labels, ...)

    ## slope of arrows
    b <- X[,2] / X[,1]

    ## offset based on string dimensions
    off <- cbind(sign(X[,1]) * (w/2 + h/4), 0.75 * h * sign(X[,2]))

    ## move the centre of the string to the continuation of the arrow
    for(i in seq_len(nrow(X))) {
        move <- off[i,2] / b[i]
        ## arrow points to the top/bottom of the text box
        if (is.finite(move) && abs(move) <= abs(off[i, 1]))
            off[i, 1] <- move
        else {
            ## arrow points to a side of the text box
            move <- b[i] * off[i,1]
            off[i, 2] <- move
        }
    }
    off + X
}
