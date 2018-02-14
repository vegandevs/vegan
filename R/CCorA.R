`CCorA` <-
    function(Y, X, stand.Y = FALSE, stand.X = FALSE, permutations = 0, ...)
{
    epsilon <- sqrt(.Machine$double.eps)
    ##
    ## BEGIN: Internal functions
    ##
    cov.inv <- function(mat, no, epsilon) {
        ## This function returns:
        ## 1) mat = matrix F of the principal components (PCA object scores);
        ## 2) S.inv = the inverse of the covariance matrix;
        ## 3) m = the rank of matrix 'mat'
        ## The inverse of the PCA covariance matrix is the diagonal
        ## matrix of (1/eigenvalues). If ncol(mat) = 1, the
        ## inverse of the covariance matrix contains 1/var(mat).
        mat <- as.matrix(mat)    # 'mat' was centred before input to cov.inv
        if(ncol(mat) == 1) {
            S.inv <- as.matrix(1/var(mat))
            m <- 1
        } else {
            S.svd <- svd(cov(mat))
            m <- ncol(mat)
            mm <- length(which(S.svd$d > max(epsilon, epsilon * S.svd$d[1L])))
            if(mm < m) {
                message(gettextf("matrix %d: rank=%d < order %d",
                                 no, mm, m))
                m <- mm
            }
            S.inv <- diag(1/S.svd$d[1:m])
            mat <- mat %*% S.svd$u[,1:m] # S.svd$u = normalized eigenvectors
        }
        list(mat=mat, S.inv=S.inv, m=m)
    }
    ## Check zero variances
    var.null <- function (mat, no) {
        problems <- diag(cov(mat)) <= 0
        if (any(problems)) {
            whichProbs <- paste(which(problems), collapse=", ")
            warning("zero variance in variable(s) ", whichProbs)
            stop("verify/modify your matrix No. ", no)
        }
        invisible(0)
    }
    probPillai <- function(Y.per, X, n, S11.inv, S22.inv, s, df1, df2, epsilon,
                           Fref, permat, ...) {
        ## Permutation test for Pillai's trace in CCorA.
        ## Reference: Brian McArdle's unpublished graduate course notes.
        S12.per <- cov(Y.per,X)
        gross.mat <- S12.per %*% S22.inv %*% t(S12.per) %*% S11.inv
        Pillai.per <- sum(diag(gross.mat))
        Fper  <- (Pillai.per*df2)/((s-Pillai.per)*df1)
        Fper >= (Fref-epsilon)
    }
    ## END: internal functions
    ##
    Y <- as.matrix(Y)
    var.null(Y,1)
    nY <- nrow(Y)
    p <- ncol(Y)
    if(is.null(colnames(Y))) {
        Ynoms <- paste("VarY", 1:p, sep="")
        } else {
        Ynoms <- colnames(Y)
        }
    X <- as.matrix(X)
    var.null(X,2)
    nX <- nrow(X)
    q <- ncol(X)
    if(is.null(colnames(X))) {
        Xnoms <- paste("VarX", 1:q, sep="")
        } else {
        Xnoms <- colnames(X)
        }
    if(nY != nX) stop("different numbers of rows in Y and X")
    n <- nY
    if(is.null(rownames(X)) & is.null(rownames(Y))) {
        rownoms <- paste("Obj", 1:n, sep="")
        } else {
        if(is.null(rownames(X))) {
            rownoms <- rownames(Y)
            } else {
            rownoms <- rownames(X)
            }
        }
    Y.c <- scale(Y, center = TRUE, scale = stand.Y)
    X.c <- scale(X, center = TRUE, scale = stand.X)
    ## Check for identical matrices
    if(p == q) {
    	if(sum(abs(Y-X)) < epsilon^2) stop("Y and X are identical")
    	if(sum(abs(Y.c-X.c)) < epsilon^2) stop("after centering, Y and X are identical")
    	}
    ## Replace Y.c and X.c by tables of their PCA object scores, computed by SVD
    temp <- cov.inv(Y.c, 1, epsilon)
    Y <- temp$mat
    pp <- temp$m
    rownames(Y) <- rownoms
    temp <- cov.inv(X.c, 2, epsilon)
    X <- temp$mat
    qq <- temp$m
    rownames(X) <- rownoms
    ## Correction PL, 26dec10
    if(max(pp,qq) >= (n-1))
    	stop("not enough degrees of freedom: max(pp,qq) >= (n-1)")
    ## Covariance matrices, etc. from the PCA scores
    S11 <- cov(Y)
    if(sum(abs(S11)) < epsilon) return(0)
    S22 <- cov(X)
    if(sum(abs(S22)) < epsilon) return(0)
    S12 <- cov(Y,X)
    if(sum(abs(S12)) < epsilon) return(0)
    S11.chol <- chol(S11)
    S11.chol.inv <- solve(S11.chol)
    S22.chol <- chol(S22)
    S22.chol.inv <- solve(S22.chol)
    ## K summarizes the correlation structure between the two sets of variables
    K <- t(S11.chol.inv) %*% S12 %*% S22.chol.inv
    K.svd <- svd(K)
    Eigenvalues <- K.svd$d^2
    ##
    ## Check for circular covariance matrix
    if((p == q) & (var(K.svd$d) < epsilon))
    	warning("[nearly] circular covariance matrix - the solution may be meaningless")
    ## K.svd$u %*% diag(K.svd$d) %*% t(K.svd$v)   # To check that K = U D V'
    axenames <- paste("CanAxis",seq_along(K.svd$d),sep="")
    U <- K.svd$u
    V <- K.svd$v
    A <- S11.chol.inv %*% U
    B <- S22.chol.inv %*% V
    Cy <- (Y %*% A)    # Correction 27dec10: remove /sqrt(n-1)
    Cx <- (X %*% B)    # Correction 27dec10: remove /sqrt(n-1)
    ## Compute the 'Biplot scores of Y and X variables' a posteriori --
    corr.Y.Cy <- cor(Y.c, Cy)  # To plot Y in biplot in space Y
    corr.Y.Cx <- cor(Y.c, Cx)  # Available for plotting Y in space of X
    corr.X.Cy <- cor(X.c, Cy)  # Available for plotting X in space of Y
    corr.X.Cx <- cor(X.c, Cx)  # To plot X in biplot in space X
    ## Add row and column names
    rownames(Cy) <- rownames(Cx) <- rownoms
	colnames(Cy) <- colnames(Cx) <- axenames
    rownames(corr.Y.Cy) <- rownames(corr.Y.Cx) <- Ynoms
    rownames(corr.X.Cy) <- rownames(corr.X.Cx) <- Xnoms
    colnames(corr.Y.Cy) <- colnames(corr.Y.Cx) <- axenames
    colnames(corr.X.Cy) <- colnames(corr.X.Cx) <- axenames

    ## Compute the two redundancy statistics
    RsquareY.X <- simpleRDA2(Y, X)
    RsquareX.Y <- simpleRDA2(X, Y)
    Rsquare.adj.Y.X <- RsquareAdj(RsquareY.X$Rsquare, n, RsquareY.X$m)
    Rsquare.adj.X.Y <- RsquareAdj(RsquareX.Y$Rsquare, n, RsquareX.Y$m)
    ## Compute Pillai's trace = sum of the canonical eigenvalues
    ##                        = sum of the squared canonical correlations
    S11.inv <- S11.chol.inv %*% t(S11.chol.inv)
    S22.inv <- S22.chol.inv %*% t(S22.chol.inv)
    gross.mat <- S12 %*% S22.inv %*% t(S12) %*% S11.inv
    PillaiTrace <- sum(diag(gross.mat))
    s   <- min(pp, qq)
    df1 <- max(pp,qq)
    df2 <- (n - max(pp,qq) - 1)
    Fval  <- (PillaiTrace*df2)/((s-PillaiTrace)*df1)
    p.Pillai <- pf(Fval, s*df1, s*df2, lower.tail=FALSE)
    permat <- getPermuteMatrix(permutations, n, ...)
    nperm <- nrow(permat)
    if (ncol(permat) != n)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                      ncol(permat), n))

    if (nperm > 0) {
        p.perm <- sapply(seq_len(nperm), function(indx, ...)
                         probPillai(Y[permat[indx,],] , X, n, S11.inv, S22.inv, s,
                                    df1, df2, epsilon, Fval, nperm, ...))
        p.perm <- (sum(p.perm) +1)/(nperm + 1)
    } else {
        p.perm <- NA
    }

    out <- list(Pillai=PillaiTrace, Eigenvalues=Eigenvalues, CanCorr=K.svd$d,
                Mat.ranks=c(RsquareX.Y$m, RsquareY.X$m),
                RDA.Rsquares=c(RsquareY.X$Rsquare, RsquareX.Y$Rsquare),
                RDA.adj.Rsq=c(Rsquare.adj.Y.X, Rsquare.adj.X.Y),
                nperm=nperm, p.Pillai=p.Pillai, p.perm=p.perm, Cy=Cy, Cx=Cx,
                corr.Y.Cy=corr.Y.Cy, corr.X.Cx=corr.X.Cx, corr.Y.Cx=corr.Y.Cx,
                corr.X.Cy=corr.X.Cy, control = attr(permat, "control"),
                call = match.call())
    class(out) <- "CCorA"
    out
}
