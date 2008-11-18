## S3 summary method for permat
`summary.permat` <-
function(object, ...)
{
    x <- object
    n <- attr(x, "times")
    if (attr(x, "ptype") != "sar" && !is.null(x$specs$reg) || !is.null(x$specs$hab))
        restr <- TRUE else restr <- FALSE
    if (restr) {
        if (!is.null(x$specs$reg) && is.null(x$specs$hab)) int <- x$specs$reg
        if (is.null(x$specs$reg) && !is.null(x$specs$hab)) int <- x$specs$hab
        if (!is.null(x$specs$reg) && !is.null(x$specs$hab))
            int <- interaction(x$specs$reg, x$specs$hab, drop=TRUE)
	nlev <- length(unique(int))        
	ssum <- numeric(n)}
    bray <- psum <- pfill <- vrow <- vcol <- numeric(n)
    for (i in 1:n) {
        bray[i] <- sum(abs(x$orig-x$perm[[i]]))/sum(x$orig+x$perm[[i]])
        psum[i] <- sum(x$orig) == sum(x$perm[[i]])
        pfill[i] <- sum(x$orig > 0) == sum(x$perm[[i]] > 0)
        vrow[i] <- sum(rowSums(x$orig) == rowSums(x$perm[[i]])) == nrow(x$orig)
        vcol[i] <- sum(colSums(x$orig) == colSums(x$perm[[i]])) == ncol(x$orig)
        if (restr) ssum[i] <- {sum(rowSums(aggregate(x$orig,list(int),sum)[,-1]) ==
            rowSums(aggregate(x$perm[[i]],list(int),sum)[,-1])) == nlev}
        }
    strsum <- if (restr) sum(ssum)/n else NA
    test <- c(sum=sum(psum)/n, fill=sum(pfill)/n, rowsums=sum(vrow)/n, colsums=sum(vcol)/n, strsum=strsum)
    x$perm <- NULL
    out <- list(x=x, bray=bray, test=test, restr=restr)
    class(out) <- c("summary.permat", "list")
    return(out)
}
