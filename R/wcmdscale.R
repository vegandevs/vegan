`wcmdscale` <-
function(d, k, eig = FALSE, add = FALSE, x.ret = FALSE, w)
{
    weight.centre <- function(x, w) {
        w.c <- apply(x, 2, weighted.mean, w = w)
        x <- sweep(x, 2, w.c, "-")
        x
    }
    if (add)
        .NotYetUsed("add")
    ZERO <- sqrt(.Machine$double.eps)
    if (!inherits(d, "dist")) {
        op <- options(warn = 2)
        on.exit(options(op))
        d <- as.dist(d)
        options(op)
    }
    m <- as.matrix(d^2)
    n <- nrow(m)
    if (missing(w))
        w <- rep(1, n)
    m <- weight.centre(m, w)
    m <- t(weight.centre(t(m), w))
    m <- m * sqrt(w) %o% sqrt(w)
    e <- eigen(-m/2, symmetric = TRUE)
    if (missing(k))
        k <- sum(e$values > ZERO)
    ev <- e$values[1:k]
    points <- sweep(e$vectors[, 1:k, drop=FALSE], 1, sqrt(w), "/")
    points <- sweep(points, 2, sqrt(ev), "*")
    rownames(points) <- rownames(m)
    if (eig || x.ret || add) {
        out <- list(points = points, eig = if (eig) e$values[-n],
                    x = if (x.ret) m, ac = NA, GOF = NA, weights = w)
        class(out) <- "wcmdscale"
    }
    else out <- points
    out
}
