## S3 summary method for permat
`summary.permat` <-
function(object, ...)
{
    x <- object
    n <- attr(x, "times")
    ss <- sum(x$orig)
    fi <- sum(x$orig > 0)
    rs <- rowSums(x$orig)
    cs <- colSums(x$orig)
    nr <- nrow(x$orig)
    nc <- ncol(x$orig)
    bray <- sapply(x$perm, function(z) sum(abs(x$orig - z)) / sum(x$orig + z))
    psum <- sapply(x$perm, function(z) ss == sum(z))
    pfill <- sapply(x$perm, function(z) fi == sum(z > 0))
    vrow <- sapply(x$perm, function(z) sum(rs == rowSums(z)) == nr)
    vcol <- sapply(x$perm, function(z) sum(cs == colSums(z)) == nc)
    if (attr(x, "is.strat")) {
        int <- attr(x, "strata")
        nlev <- length(unique(int))
        rsagg <- rowSums(aggregate(x$orig, list(int), sum)[,-1])
        ssum <- sapply(x$perm, function(z)
            sum(rsagg == rowSums(aggregate(z, list(int), sum)[,-1])) == nlev)
    } else ssum <- NULL
    x$perm <- NULL
    out <- list(x=x, bray=bray, sum=psum, fill=pfill, rowsums=vrow, colsums=vcol, strsum=ssum)
    class(out) <- c("summary.permat", "list")
    return(out)
}
