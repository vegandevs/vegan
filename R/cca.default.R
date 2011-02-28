`cca.default` <-
    function (X, Y, Z, ...) 
{
    ZERO <- 1e-04
    CCA <- NULL
    pCCA <- NULL
    CA <- NULL
    weight.centre <- function(x, w) {
        w.c <- apply(x, 2, weighted.mean, w = w)
        x <- sweep(x, 2, w.c, "-")
        x <- sweep(x, 1, sqrt(w), "*")
        attr(x, "centre") <- w.c
        x
    }
    X <- as.matrix(X)
    if (any(rowSums(X) <= 0)) 
        stop("All row sums must be >0 in the community data matrix")
    if (any(tmp <- colSums(X) <= 0)) {
        exclude.spec <- seq(along=tmp)[tmp]
        names(exclude.spec) <- colnames(X)[tmp]
        class(exclude.spec) <- "exclude"
        X <- X[, !tmp, drop = FALSE]
    }
    gran.tot <- sum(X)
    X <- X/gran.tot
    rowsum <- apply(X, 1, sum)
    colsum <- apply(X, 2, sum)
    rc <- outer(rowsum, colsum)
    Xbar <- (X - rc)/sqrt(rc)
    tot.chi <- sum(svd(Xbar, nu = 0, nv = 0)$d^2)
    if (!missing(Z) && !is.null(Z)) {
        Z <- as.matrix(Z)
        Z.r <- weight.centre(Z, rowsum)
        Q <- qr(Z.r)
        Z <- qr.fitted(Q, Xbar)
        tmp <- sum(svd(Z, nu = 0, nv = 0)$d^2)
        if (Q$rank) {
            pCCA <- list(rank = Q$rank, tot.chi = tmp, QR = Q, 
                         Fit = Z, envcentre = attr(Z.r, "centre"))
            Xbar <- qr.resid(Q, Xbar)
        }
        if (tmp < ZERO)
            pCCA$tot.chi <- 0
    }
    else Z.r <- NULL
    if (!missing(Y) && !is.null(Y)) {
        Y <- as.matrix(Y)
        Y.r <- weight.centre(Y, rowsum)
        Q <- qr(cbind(Z.r, Y.r), tol = ZERO)
        if (is.null(pCCA)) 
            rank <- Q$rank
        else rank <- Q$rank - pCCA$rank
        ## save rank of constraints
        qrank <- rank
        Y <- qr.fitted(Q, Xbar)
        sol <- svd(Y)
        ## rank of svd can be < qrank
        rank <- min(rank, sum(sol$d > ZERO))
        ax.names <- paste("CCA", 1:length(sol$d), sep = "")
        colnames(sol$u) <- ax.names
        colnames(sol$v) <- ax.names
        names(sol$d) <- ax.names
        rownames(sol$u) <- rownames(X)
        rownames(sol$v) <- colnames(X)
        if (rank) {
            CCA <- list(eig = sol$d[1:rank]^2)
            CCA$u <- sweep(as.matrix(sol$u[, 1:rank, drop = FALSE]), 
                           1, 1/sqrt(rowsum), "*")
            CCA$v <- sweep(as.matrix(sol$v[, 1:rank, drop = FALSE]), 
                           1, 1/sqrt(colsum), "*")
            CCA$u.eig <- sweep(CCA$u, 2, sol$d[1:rank], "*")
            CCA$v.eig <- sweep(CCA$v, 2, sol$d[1:rank], "*")
            CCA$wa.eig <- sweep(Xbar %*% sol$v[, 1:rank, drop = FALSE], 
                                1, 1/sqrt(rowsum), "*")
            CCA$wa <- sweep(CCA$wa.eig, 2, 1/sol$d[1:rank], "*")
            oo <- Q$pivot
            if (!is.null(pCCA$rank)) 
                oo <- oo[-(1:pCCA$rank)] - ncol(Z.r)
            oo <- oo[1:qrank]
            if (length(oo) < ncol(Y.r)) 
                CCA$alias <- colnames(Y.r)[-oo]
            CCA$biplot <- cor(Y.r[, oo, drop = FALSE], sol$u[, 
                                        1:rank, drop = FALSE])
            CCA$rank <- rank
            CCA$qrank <- qrank
            CCA$tot.chi <- sum(CCA$eig)
            CCA$QR <- Q
            CCA$envcentre <- attr(Y.r, "centre")
            CCA$Xbar <- Xbar
        } else {                # zero rank
            CCA <- list(eig = 0, rank = rank, qrank = qrank, tot.chi = 0,
                        QR = Q, Xbar = Xbar)
            u <- matrix(0, nrow=nrow(sol$u), ncol=0)
            v <- matrix(0, nrow=nrow(sol$v), ncol=0)
            CCA$u <- CCA$u.eig <- CCA$wa <- CCA$wa.eig <- u
            CCA$v <- CCA$v.eig <- v
            CCA$biplot <- matrix(0, 0, 0)
            CCA$alias <- colnames(Y.r)
        }
        Xbar <- qr.resid(Q, Xbar)
        if (exists("exclude.spec")) {
            attr(CCA$v, "na.action") <- exclude.spec
            attr(CCA$v.eig, "na.action") <- exclude.spec
        }
        
    }
    Q <- qr(Xbar)
    sol <- svd(Xbar)
    ax.names <- paste("CA", 1:length(sol$d), sep = "")
    colnames(sol$u) <- ax.names
    colnames(sol$v) <- ax.names
    names(sol$d) <- ax.names
    rownames(sol$u) <- rownames(X)
    rownames(sol$v) <- colnames(X)
    rank <- min(Q$rank, sum(sol$d > ZERO))
    if (rank) {
        CA <- list(eig = sol$d[1:rank]^2)
        CA$u <- sweep(as.matrix(sol$u[, 1:rank, drop = FALSE]), 
                      1, 1/sqrt(rowsum), "*")
        CA$v <- sweep(as.matrix(sol$v[, 1:rank, drop = FALSE]), 
                      1, 1/sqrt(colsum), "*")
        CA$u.eig <- sweep(CA$u, 2, sol$d[1:rank], "*")
        CA$v.eig <- sweep(CA$v, 2, sol$d[1:rank], "*")
        CA$rank <- rank
        CA$tot.chi <- sum(CA$eig)
        CA$Xbar <- Xbar
        
    } else {   # zero rank: no residual component
        CA <- list(eig = 0, rank = rank, tot.chi = 0,
                   Xbar = Xbar)
        CA$u <- CA$u.eig <- matrix(0, nrow(sol$u), 0)
        CA$v <- CA$v.eig <- matrix(0, nrow(sol$v), 0)
    }
    if (exists("exclude.spec")) {
        attr(CA$v, "na.action") <- exclude.spec
        attr(CA$v.eig, "na.action") <- exclude.spec
    }
    call <- match.call()
    call[[1]] <- as.name("cca")
    ## computed pCCA$rank was needed before, but zero it here
    if (!is.null(pCCA) && pCCA$tot.chi == 0)
        pCCA$rank <- 0
    sol <- list(call = call, grand.total = gran.tot, rowsum = rowsum, 
                colsum = colsum, tot.chi = tot.chi, pCCA = pCCA, CCA = CCA, 
                CA = CA)
    sol$method <- "cca"
    sol$inertia <- "mean squared contingency coefficient"
    class(sol) <- "cca"
    sol
}
