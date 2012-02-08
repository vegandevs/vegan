indpower <-
function(x, type=0)
{
    x <- as.matrix(x)
    x <- ifelse(x > 0, 1, 0)
    if (NCOL(x) < 2)
        stop("provide at least 2 columns for 'x'")
    if (!(type %in% 0:2))
        stop("'type' must be in c(0, 1, 2)")
    n <- nrow(x)
    j <- t(x) %*% x
    ip1 <- sweep(j, 1, diag(j), "/")
    ip2 <- 1 - sweep(-sweep(j, 2, diag(j), "-"), 1, n - diag(j), "/")
    out <- switch(as.character(type),
        "0" = sqrt(ip1 * ip2),
        "1" = ip1,
        "2" = ip2)
    cn <- if (is.null(colnames(out)))
        1:ncol(out) else colnames(out)
    rn <- if (is.null(rownames(out)))
        1:ncol(out) else rownames(out)
    colnames(out) <- paste("t", cn, sep=".")
    rownames(out) <- paste("i", rn, sep=".")
    out
}
