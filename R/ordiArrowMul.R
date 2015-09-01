### Scaling of arrows to 'fill' a plot with vectors centred at 'at'.
### Plot dims from 'par("usr")' and arrow heads are in 'x'.
`ordiArrowMul` <- function (x, at = c(0,0), fill = 0.75,
                            display, choices = c(1,2), ...) {
    ## handle x, which we try with scores, but also retain past usage of
    ## a two column matrix
    X <- if (is.matrix(x)) {
        nc <- NCOL(x)
        if (nc != 2L) {
            stop("A 2-column matrix of coordinates is required & not supplied.")
        }
        x
    } else {
        if (inherits(x, "envfit")) {
            scores(x, display = "vectors", ...)[, 1:2]
        } else {
            scores(x, display = display, choices = choices, ...)
        }
    }

    u <- par("usr")
    u <- u - rep(at, each = 2)
    r <- c(range(X[,1], na.rm = TRUE), range(X[,2], na.rm = TRUE))
    ## 'rev' takes care of reversed axes like xlim(1,-1)
    rev <- sign(diff(u))[-2]
    if (rev[1] < 0)
        u[1:2] <- u[2:1]
    if (rev[2] < 0)
        u[3:4] <- u[4:3]
    u <- u/r
    u <- u[is.finite(u) & u > 0]
    fill * min(u)
}
