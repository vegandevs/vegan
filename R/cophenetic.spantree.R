`cophenetic.spantree` <-
    function(x)
{
    n <- x$n
    mat <- matrix(NA, nrow=n, ncol=n)
    if (n < 2)
        return(as.dist(mat))
    ind <- apply(cbind(2:n, x$kid), 1, sort)
    ind <- t(ind[2:1,])
    mat[ind] <- x$dist
    d <- as.dist(mat)
    attr(d, "Labels") <- x$labels
    stepacross(d, path = "extended", toolong=0, trace=FALSE)
}
