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
    ## Remove zero eigenvalues, keep negative
    keep <- abs(e$values) > ZERO
    e$values <- e$values[keep]
    e$vectors <- e$vectors[, keep, drop = FALSE]
    ## Deweight and scale axes -- also negative
    points <- sweep(e$vectors, 1, sqrt(w), "/")
    points <- sweep(points, 2, sqrt(abs(e$values)), "*")
    rownames(points) <- rownames(m)
    ## If 'k' not given, find it as the number of positive
    ## eigenvalues, and also save negative eigenvalues
    negaxes <- NULL
    if (missing(k) || k > sum(e$value > ZERO)) {
        k <- sum(e$values > ZERO)
        if (any(e$values < 0))
            negaxes <- points[, e$values < 0, drop = FALSE]
    }
    points <- points[, 1:k, drop=FALSE]
    if (eig || x.ret || add) {
        out <- list(points = points, eig = if (eig) e$values,
                    x = if (x.ret) m, ac = NA, GOF = NA, weights = w,
                    negaxes = negaxes)
        class(out) <- "wcmdscale"
    }
    else out <- points
    out
}
