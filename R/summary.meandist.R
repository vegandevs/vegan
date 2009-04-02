`summary.meandist` <-
    function(object, ...)
{
    n <- attr(object, "n")
    wmat <- n %o% n
    diag(wmat) <- diag(wmat) - n
    ## mean distances within, between groups and in total
    W <- weighted.mean(diag(object), w = diag(wmat), na.rm = TRUE)
    B <- weighted.mean(object[lower.tri(object)],
                       w = wmat[lower.tri(wmat)], na.rm = TRUE)
    D <- weighted.mean(object, w = wmat, na.rm = TRUE)
    ## Variants of MRPP statistics
    A1 <- weighted.mean(diag(object), w = n, na.rm = TRUE)
    A2 <- weighted.mean(diag(object), w = n - 1, na.rm = TRUE)
    A3 <- weighted.mean(diag(object), w = n * (n - 1), na.rm = TRUE)
    ##
    out <- list(W = W, B = B, D = D, CS = B-W,
                A1 = 1 - A1/D, A2 = 1 - A2/D, A3 = 1 - A3/D)
    class(out) <- "summary.meandist"
    out
}

