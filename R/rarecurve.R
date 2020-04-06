`rarecurve` <-
    function(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
             label = TRUE, col, lty, ...)
{
    ## matrix is faster than data.frame
    x <- as.matrix(x)
    ## check input data: must be counts
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    ## sort out col and lty
    if (missing(col))
        col <- par("col")
    if (missing(lty))
        lty <- par("lty")
    tot <- rowSums(x)
    S <- specnumber(x)
    ## remove empty rows or we fail
    if (any(S <= 0)) {
        message("empty rows removed")
        x <- x[S > 0,, drop =FALSE]
        tot <- tot[S > 0]
        S <- S[S > 0]
    }
    nr <- nrow(x)
    ## rep col and lty to appropriate length
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
    ## Rarefy
    out <- lapply(seq_len(nr), function(i) {
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
            ## don't want names on n an `c` adds a name from `tot[i]`)
            n <- c(n, tot[i], use.names = FALSE)
        }
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
    for (ln in seq_along(out)) {
        N <- attr(out[[ln]], "Subsample")
        lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
    }
    ## label curves at their endpoitns
    if (label) {
        ordilabel(cbind(tot, S), labels=rownames(x), ...)
    }
    invisible(out)
}
