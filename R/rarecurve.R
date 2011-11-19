`rarecurve` <-
    function(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
             label = TRUE,...)
{
    tot <- rowSums(x)
    S <- specnumber(x)
    nr <- nrow(x)
    ## Rarefy
    out <- lapply(seq_len(nr), function(i) {
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i])
            n <- c(n, tot[i])
        drop(rarefy(x[i,], n))
    })
    Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
    Smax <- sapply(out, max)
    ## set up plot
    plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
         type = "n", ...)
    ## rarefied richnesses for given 'sample'
    if (!missing(sample)) {
        abline(v = sample)
        rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), y = z,
                                             xout = sample, rule = 1)$y)
        abline(h = rare, lwd=0.5)
    }
    ## rarefaction curves
    for(ln in seq_len(length(out))) {
        N <- attr(out[[ln]], "Subsample")
        lines(N, out[[ln]], ...)
    }
    ## label curves at their endpoitns
    if (label) {
        ordilabel(cbind(tot, S), labels=rownames(x), ...)
    }
    invisible(out)
}
