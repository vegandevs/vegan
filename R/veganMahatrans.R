### Internal function for Mahalanobis transformation of the matrix.
### Mahalanobis transformation of matrix X is M = X S^(-1/2) where S
### is the covariance matrix. The inverse square root of S is found
### via eigen decomposition S = G L G^T, where G is the matrix of
### eigenvectors, and L is the diagonal matrix of eigenvalues. Thus
### S^(-1/2) = G L^(-1/2) G^T. This is an internal function so that
### input must be correct: 'x' must be a centred matrix (not a
### data.frame, not raw data).
`veganMahatrans` <-
    function (x, s2, tol = sqrt(.Machine$double.eps), na.rm = FALSE)
{
    if (missing(s2))
        s2 <- cov(x, use = if(na.rm) "pairwise.complete.obs" else "all.obs")
    e <- eigen(s2, symmetric = TRUE)
    k <- e$values > max(tol, tol * e$values[1L])
    sisqr <- e$vectors[,k, drop=FALSE] %*%
        (sqrt(1/e$values[k]) * t(e$vectors[,k, drop = FALSE]))
    x %*% sisqr
}
