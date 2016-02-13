`wascores` <-
    function (x, w, expand = FALSE) 
{
    if(any(w < 0) || sum(w) == 0)
        stop("weights must be non-negative and not all zero")
    x <- as.matrix(x)
    w <- as.matrix(w)
    nc <- ncol(x)
    nr <- ncol(w)
    wa <- matrix(NA, nrow = nr, ncol = nc)
    colnames(wa) <- colnames(x)
    rownames(wa) <- colnames(w)
    for (i in 1:nr) {
        wa[i, ] <- apply(x, 2, weighted.mean, w = w[, i])
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
        attr(wa, "shrinkage") <- 1/mul^2
        attr(wa, "centre") <- wa.cov$center
    }
    wa
}
