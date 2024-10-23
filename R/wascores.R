`wascores` <-
    function (x, w, expand = FALSE, stdev = FALSE)
{
    if(any(w < 0) || sum(w) == 0)
        stop("weights must be non-negative and not all zero")
    x <- as.matrix(x)
    w <- as.matrix(w)
    nc <- ncol(x)
    nr <- ncol(w)
    wa <- t(sapply(colnames(w),
                   function(m) apply(x, 2, weighted.mean, w = w[,m])))
    if (stdev) {
        sdwa <- sqrt(sapply(colnames(x), function(k)
            sapply(colnames(w), function(m)
                cov.wt(x[,k,drop=FALSE], w[,m])$cov)))
    }
    if (expand) {
        i <- complete.cases(wa)
        x.w <- rowSums(w)
        ewa.w <- colSums(w[,i, drop=FALSE])
        ewa <- wa[i,, drop=FALSE]
        x.cov <- cov.wt(x, x.w, method = "ML")
        wa.cov <- cov.wt(ewa, ewa.w, method = "ML")
        mul <- sqrt(diag(x.cov$cov)/diag(wa.cov$cov))
        ewa <- sweep(ewa, 2, wa.cov$center, "-")
        ewa <- sweep(ewa, 2, mul, "*")
        ewa <- sweep(ewa, 2, wa.cov$center, "+")
        wa[i,] <- ewa
        if (stdev)
            sdwa <- sweep(sdwa, 2, mul, "*")
        attr(wa, "shrinkage") <- 1/mul^2
        attr(wa, "centre") <- wa.cov$center
    }
    if (stdev) {
        wa <- list("wa" = wa, "stdev" = sdwa,
                   "N2" = diversity(w, "invsimpson", MARGIN=2))
    }
    wa
}
