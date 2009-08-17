`treedist` <-
    function(x, tree, ...)
{
    n <- nrow(x)
    ABJ <- matrix(0, n , n)
    dmat <- as.matrix(cophenetic(tree))
    for(j in 1:n)
        for (k in j:n) {
            jk <- x[j,] > 0 | x[k,] > 0
            if (sum(jk) > 1)
                ABJ[k, j] <- treeheight(update(tree, d = as.dist(dmat[jk, jk]))) 
        }
    A <- diag(ABJ)
    AB <- as.dist(outer(A, A, "+"))
    ABJ <- as.dist(ABJ)
    out <- (2 * ABJ - AB)/ABJ
    attr(out, "method") <- "treedist"
    attr(out, "call") <- match.call()
    attr(out, "Labels") <- row.names(x)
    out
}
