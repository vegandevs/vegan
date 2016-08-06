`rarefy` <-
    function (x, sample, se = FALSE, MARGIN = 1) 
{
    x <- as.matrix(x)
    ## as.matrix changes an n-vector to a n x 1 matrix
    if (ncol(x) == 1 && MARGIN == 1)
        x <- t(x)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    minsample <- min(apply(x, MARGIN, sum))
    if (missing(sample)) {
        stop(
            gettextf(
                "The size of 'sample' must be given --\nHint: Smallest site maximum %d",
                minsample))
    }
    if (any(sample > minsample))
        warning(
            gettextf(
                "Requested 'sample' was larger than smallest site maximum (%d)",
                minsample))
    rarefun <- function(x, sample) {
        x <- x[x > 0]
        J <- sum(x)
        ldiv <- lchoose(J, sample)
        p1 <- ifelse(J - x < sample, 0, exp(lchoose(J - x, sample) - 
                                            ldiv))
        out <- sum(1 - p1)
        if (se) {
            V <- sum(p1 * (1 - p1))
            Jxx <- J - outer(x, x, "+")
            ind <- lower.tri(Jxx)
            Jxx <- Jxx[ind]
            V <- V + 2 * sum(ifelse(Jxx < sample, 0, exp(lchoose(Jxx, 
                                                                 sample) - ldiv)) - outer(p1, p1)[ind])
            ## V is >= 0, but numerical zero can be negative (e.g,
            ## -1e-16), and we avoid taking its square root
            out <- cbind(out, sqrt(max(V, 0)))
        }
        out
    }
    if (length(sample) > 1) {
        S.rare <- sapply(sample, function(n) apply(x, MARGIN, rarefun, sample = n))
        S.rare <- matrix(S.rare, ncol=length(sample))
        colnames(S.rare) <- paste("N", sample, sep="")
        if (se) {
            dn <- unlist(dimnames(x)[MARGIN])
            rownames(S.rare) <- paste(rep(dn, each=2), c("S","se"), sep=".")
        }
    } else {
        S.rare <- apply(x, MARGIN, rarefun, sample = sample)
        if (se) 
            rownames(S.rare) <- c("S", "se")
    }
    attr(S.rare, "Subsample") <- sample
    S.rare
}
