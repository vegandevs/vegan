`CCorA` <-
    function(Y, X, stand.Y = FALSE, stand.X = FALSE, nperm = 0, ...)
{
    require(MASS) || stop("requires packages 'MASS'")
    epsilon <- sqrt(.Machine$double.eps)
    ## BEGIN: Internal functions
    cov.inv <- function(mat, no, epsilon) {
        ## This function returns:
        ## 1) mat = the matrix of PCA object scores (by SVD);
        ## 2) S.inv = the inverse of the covariance matrix;
        ## 3) m = the rank of matrix 'mat'
        ## The inverse of the PCA covariance matrix is simply the
        ## diagonal matrix of (1/eigenvalues).  If ncol(mat) = 1, the
        ## inverse of the covariance matrix simply contains
        ## 1/var(mat).
        mat <- as.matrix(mat)
        if(ncol(mat) == 1) {
            S.inv <- as.matrix(1/var(mat))
            m <- 1
        } else {
            S.svd <- svd(cov(mat))
            m <- ncol(mat)
            mm <- length(which(S.svd$d > epsilon))
            if(mm < m) {
                message("Information - Matrix",no,": rank=",mm," < order",m)
                m <- mm
            }
            S.inv <- diag(1/S.svd$d[1:m])
            mat <- mat %*% S.svd$u[,1:m]
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
    probPillai <- function(Y, X, n, S11.inv, S22.inv, s, df1, df2, epsilon,
                           Fref, nperm, ...) {
        ## Permutation test for Pillai's trace in CCorA.
        ## Reference: Brian McArdle's unpublished graduate course notes.
        nGE <- 1
        for(i in 1:nperm) {
            Y.per <- Y[permuted.index(n, ...),, drop=FALSE]
            S12.per <- cov(Y.per,X)
            gross.mat <- S12.per %*% S22.inv %*% t(S12.per) %*% S11.inv
            Pillai.per <- sum(diag(gross.mat))
            Fper  <- (Pillai.per*df2)/((s-Pillai.per)*df1)
            if(Fper >= (Fref-epsilon)) nGE <- nGE+1
        }
        P <- nGE/(nperm+1)
    }
    ## END: internal functions
    Y <- as.matrix(Y)
    var.null(Y,1)
    nY <- nrow(Y)
    p <- ncol(Y)
    Ynoms <- colnames(Y)
    X <- as.matrix(X)
    var.null(X,2)
    nX <- nrow(X)
    q <- ncol(X)
    Xnoms <- colnames(X)
    if(nY != nX) stop("Different numbers of rows in Y and X")
    n <- nY
    if((p+q) >= (n-1)) stop("Not enough degrees of freedom!")
    rownoms <- rownames(Y)
    Y.c <- scale(Y, center = TRUE, scale = stand.Y)
    X.c <- scale(X, center = TRUE, scale = stand.X)
    ## Replace Y.c and X.c by tables of their PCA object scores, computed by SVD
    temp <- cov.inv(Y.c, 1, epsilon)
    Y <- temp$mat
    pp <- temp$m
    rownames(Y) <- rownoms
    temp <- cov.inv(X.c, 2, epsilon)
    X <- temp$mat
     qq <- temp$m
    rownames(X) <- rownoms
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
    EigenValues <- K.svd$d^2
    ## K.svd$u %*% diag(K.svd$d) %*% t(K.svd$v)   # To check that K = U D V'
    axenames <- paste("CanAxis",1:length(K.svd$d),sep="")
    U <- K.svd$u
    V <- K.svd$v
    A <- S11.chol.inv %*% U
    B <- S22.chol.inv %*% V
    Cy <- (Y %*% A)/sqrt(n-1)
    Cx <- (X %*% B)/sqrt(n-1)
    ## Compute the 'Biplot scores of Y variables' a posteriori --
    ## use 'ginv' for inversion in case there is collinearity
    ## AA <- coefficients of the regression of Cy on Y.c, times sqrt(n-1)
    ## AA <- sqrt(n-1) * [Y'Y]-1 Y' Cy
    YprY <- t(Y.c) %*% Y.c
    AA <- sqrt(n-1) * ginv(YprY) %*% t(Y.c) %*% Cy
    ##
    ## Compute the 'Biplot scores of X variables' a posteriori --
    XprX <- t(X) %*% X
    BB <- sqrt(n-1) * ginv(XprX) %*% t(X) %*% Cx
    ## Add row and column names
    rownames(AA) <- Ynoms
    rownames(BB) <- Xnoms
    rownames(Cy) <- rownames(Cx) <- rownoms
    colnames(AA) <- colnames(BB) <- colnames(Cy) <- colnames(Cx) <- axenames
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

    if(nperm > 0) {
        p.perm <- probPillai(Y, X, n, S11.inv, S22.inv, s, df1, df2,
                             epsilon, Fval, nperm, ...)
    } else {
        p.perm <- NA
    }

    out <- list(Pillai=PillaiTrace, EigenValues=EigenValues, CanCorr=K.svd$d,
                Mat.ranks=c(RsquareX.Y$m, RsquareY.X$m), 
                RDA.Rsquares=c(RsquareY.X$Rsquare, RsquareX.Y$Rsquare),
                RDA.adj.Rsq=c(Rsquare.adj.Y.X, Rsquare.adj.X.Y),
                nperm=nperm, p.Pillai=p.Pillai, p.perm=p.perm,
                AA=AA, BB=BB, Cy=Cy, Cx=Cx, call = match.call())
    class(out) <- "CCorA"
    out
}
