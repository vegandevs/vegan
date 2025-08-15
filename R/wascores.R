`wascores` <-
    function (x, w, expand = FALSE, stdev = FALSE)
{
    if(any(w < 0) || sum(w) == 0)
        stop("weights must be non-negative and not all zero")
    x <- as.matrix(x)
    w <- as.matrix(w)
    nc <- ncol(x)
    nr <- ncol(w)
    dnam <- list(colnames(w), colnames(x))
    wa <- t(sapply(seq_len(ncol(w)),
                   function(i) apply(x, 2, weighted.mean, w = w[,i])))
    wa <- matrix(wa, nr, nc, dimnames = dnam)
    if (stdev) {
        sdwa <- sqrt(sapply(seq_len(nc), function(k)
            sapply(colnames(w), function(m)
                cov.wt(x[,k,drop=FALSE], w[,m])$cov)))
        sdwa <- matrix(sdwa, nr, nc, dimnames = dnam)
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
        attr(wa, "call") <- match.call()
        class(wa) <- "wascores"
    }
    wa
}

`print.wascores` <-
    function(x, ...)
{
    if (!is.null(x$stdev)) { # should be always TRUE with class "wascores"
        if (ncol(x$wa) == 1) {
            out <- do.call("cbind", x)
            colnames(out) <- c("WA", "SD", "N2")
            print(out)
        } else { # more than one x-variate
            print(x$wa)
            cat("\nUse scores() to see stdev (or derived statistics) and N2\n\n")
        }
    }
    invisible(x)
}

`scores.wascores` <-
    function(x, display = c("wa", "stdev", "var", "se", "n2", "raw"), ...)
{
    display <- tolower(display)
    display <- match.arg(display)
    ## Calculation of CI via t-value is currently disabled (although
    ## there is an entry in switch). If it is ever enabled in similar
    ## way, p-value should be lifted to function arguments and "ci" in display.
    p <- 0.95
    if(display == "ci")
        tval <- qt((1 - p)/2, x$N2, lower.tail = FALSE)
    switch(display,
           "wa" = x$wa,
           "stdev" = x$stdev,
           "var" = x$stdev^2,
           "se" = x$stdev/sqrt(x$N2),
           "ci" = x$stdev/sqrt(x$N2) * tval,
           "n2" = x$N2,
           "raw" = unclass(x))
}
