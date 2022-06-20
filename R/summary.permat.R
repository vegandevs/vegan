## S3 summary method for permat
`summary.permat` <-
    function(object, ...)
{
    x <- object
    ## calculations are much faster if x$orig is matrix instead of data.frame
    x$orig <- data.matrix(x$orig)
    ss <- sum(x$orig)
    fi <- sum(x$orig > 0)
    rs <- rowSums(x$orig)
    cs <- colSums(x$orig)
    rb <- rowSums(x$orig > 0)
    cb <- colSums(x$orig > 0)
    nr <- nrow(x$orig)
    nc <- ncol(x$orig)
    bray <- sapply(x$perm, function(z) sum(abs(x$orig - z)) / sum(x$orig + z))
    psum <- sapply(x$perm, function(z) ss == sum(z))
    pfill <- sapply(x$perm, function(z) fi == sum(z > 0))
    vrow <- sapply(x$perm, function(z) sum(rs == rowSums(z)) == nr)
    vcol <- sapply(x$perm, function(z) sum(cs == colSums(z)) == nc)
    brow <- sapply(x$perm, function(z) sum(rb == rowSums(z > 0)) == nr)
    bcol <- sapply(x$perm, function(z) sum(cb == colSums(z > 0)) == nc)
    if (attr(x, "is.strat")) {
        int <- attr(x, "strata")
        nlev <- length(unique(int))
        rsagg <- rowSums(aggregate(x$orig, list(int), sum)[,-1])
        ssum <- sapply(x$perm, function(z)
            sum(rsagg == rowSums(aggregate(z, list(int), sum)[,-1])) == nlev)
    } else ssum <- NULL
    ## Chisq
    E <- rs %o% cs / ss
    chisq <- sapply(x$perm, function(z) sum((z - E)^2 / E))
    attr(chisq, "chisq.orig") <- sum((x$orig - E)^2 / E)
#    attr(chisq, "df") <- (nr - 1) * (nc - 1)
    ## ts if sequential
    seqmethods <- sapply(make.commsim(), function(z) make.commsim(z)$isSeq)
    seqmethods <- names(seqmethods)[seqmethods]
#    seqmethods <- c("swap", "tswap", "abuswap")
    if (attr(x, "method") %in% seqmethods) {
        startval <- attr(x, "burnin") + 1
        dtime <- max(1, attr(x, "thin"))
        bray <- ts(bray, start = startval, deltat = dtime)
        chisq <- ts(chisq, start = startval, deltat = dtime)
    }
    x$perm <- NULL
    out <- list(x=x, bray=bray, chisq=chisq, sum=psum, fill=pfill, rowsums=vrow, colsums=vcol,
        browsums=brow, bcolsums=bcol, strsum=ssum)
    class(out) <- c("summary.permat", "list")
    out
}
