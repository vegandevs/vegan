`rarecurve` <-
    function(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
             label = TRUE, col, lty, tidy = FALSE, ...)
{
    ## matrix is faster than data.frame
    x <- as.matrix(x, rownames.force = TRUE)
    ## check input data: must be counts
    if (!isTRUE(all.equal(x, round(x))))
        stop("function accepts only integers (counts)")
    x <- round(x) # x may not be exact integer
    ## should be observed counts
    minobs <- min(x[x > 0])
    if (minobs > 1)
        warning(gettextf("most observed count data have counts 1, but smallest count is %d", minobs))
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
        ## already warned on possibly non-observed counts: do not
        ## repeat warnings for every row
        drop(suppressWarnings(rarefy(x[i,], n)))
    })
    ## instead of plotting a rarecurve, return a "tidy" data frame and
    ## the let the user figure out how to display the results
    if (tidy) {
        len <- sapply(out, length)
        nm <- rownames(x)
        df <- data.frame(
            "Site" = factor(rep(nm, len), levels=nm),
            "Sample" = unlist(lapply(out, attr, which="Subsample")),
            "Species" = unlist(out))
        return(df) # exit with data.frame
    }
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
