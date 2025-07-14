`wcmdscale` <-
    function(d, k, eig = FALSE, add = FALSE, x.ret = FALSE, w)
{
    ## Force eig=TRUE if add, x.ret or !missing(w)
    if(x.ret)
        eig <- TRUE
    ZERO <- sqrt(.Machine$double.eps)
    if (!inherits(d, "dist")) {
        op <- options(warn = 2)
        on.exit(options(op))
        d <- as.dist(d)
        options(op)
    }
    ## handle add constant to make d Euclidean
    if (is.logical(add) && add)
        add <- "lingoes"
    if (is.character(add)) {
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- addLingoes(as.matrix(d))
            d <- sqrt(d^2 + 2 * ac)
        } else if (add == "cailliez") {
            ac <- addCailliez(as.matrix(d))
            d <- d + ac
        }
    } else {
        ac <- NA
    }
    ## Gower centring
    m <- as.matrix(d^2)
    n <- nrow(m)
    if (missing(w))
        w <- rep(1, n)
    m <- .Call(do_wcentre, m, w, PACKAGE = "vegan")
    m <- t(.Call(do_wcentre, t(m), w, PACKAGE = "vegan"))
    e <- eigen(-m/2, symmetric = TRUE)
    ## Remove zero eigenvalues, keep negative
    keep <- abs(e$values) > max(ZERO, ZERO * e$values[1L])
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
    if (k) # there may be no positive eigenvalues
        points <- points[, 1:k, drop=FALSE]
    points[!is.finite(points)] <- NA
    ## Goodness of fit
    ev <- e$values[1:k]
    ev <- ev[ev > 0]
    ## GOF for real and all axes
    GOF <- c(sum(ev)/sum(abs(e$values)),
             sum(ev)/sum(e$values[e$values > 0]))
    if (eig || x.ret) {
        if (NCOL(points) > 0)
            colnames(points) <- paste("Dim", seq_len(NCOL(points)), sep="")
        out <- list(points = points, eig = if (eig) e$values,
                    x = if (x.ret) m, ac = ac, add = add, GOF = GOF,
                    weights = w, negaxes = negaxes, call = match.call())
        class(out) <- "wcmdscale"
    }
    else out <- points
    out
}
