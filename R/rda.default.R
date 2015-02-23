`rda.default` <-
    function (X, Y, Z, scale = FALSE, ...) 
{
    ZERO <- 1e-05
    CCA <- NULL
    pCCA <- NULL
    CA <- NULL
    ## Protect against grave misuse: some people have used
    ## dissimilarities instead of data
    if (inherits(X, "dist") || NCOL(X) == NROW(X) &&
        isTRUE(all.equal(X, t(X))))
        stop("function cannot be used with (dis)similarities")
    X <- as.matrix(X)
    NR <- nrow(X) - 1
    Xbar <- scale(X, center = TRUE, scale = scale)
    SD <- apply(Xbar, 2, sd)
    if (scale) 
        Xbar[is.nan(Xbar)] <- 0
    tot.chi <- sum(svd(Xbar, nu = 0, nv = 0)$d^2)/NR
    if (!missing(Z) && !is.null(Z)) {
        Z <- as.matrix(Z)
        Z.r <- scale(Z, center = TRUE, scale = FALSE)
        Q <- qr(Z.r)
        Z <- qr.fitted(Q, Xbar)
        tmp <- sum(svd(Z, nu = 0, nv = 0)$d^2)/NR
        if (Q$rank) {
            pCCA <- list(rank = Q$rank, tot.chi = tmp, QR = Q, 
                         Fit = Z, envcentre = attr(Z.r, "scaled:center"))
            Xbar <- qr.resid(Q, Xbar)
        }
        if (tmp < ZERO)
            pCCA$tot.chi <- 0
    }
    else Z.r <- NULL
    if (!missing(Y) && !is.null(Y)) {
        Y <- as.matrix(Y)
        Y.r <- scale(Y, center = TRUE, scale = FALSE)
        Q <- qr(cbind(Z.r, Y.r), tol = ZERO)
        if (is.null(pCCA)) 
            rank <- Q$rank
        else rank <- Q$rank - pCCA$rank
        ## qrank saves the rank of the constraints
        qrank <- rank
        Y <- qr.fitted(Q, Xbar)
        sol <- svd(Y)
        ## it can happen that rank < qrank
        rank <- min(rank, sum(sol$d > (sol$d[1L] * ZERO)))
        sol$d <- sol$d/sqrt(NR)
        ax.names <- paste("RDA", seq_along(sol$d), sep = "")
        colnames(sol$u) <- ax.names
        colnames(sol$v) <- ax.names
        names(sol$d) <- ax.names
        rownames(sol$u) <- rownames(X)
        rownames(sol$v) <- colnames(X)
        if (rank) {
            CCA <- list(eig = sol$d[1:rank]^2)
            CCA$u <- as.matrix(sol$u)[, 1:rank, drop = FALSE]
            CCA$v <- as.matrix(sol$v)[, 1:rank, drop = FALSE]
            wa.eig <- Xbar %*% sol$v[, 1:rank, drop = FALSE]
            wa.eig <- wa.eig/sqrt(NR)
            CCA$wa <- sweep(wa.eig, 2, 1/sol$d[1:rank], "*")
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
            CCA$envcentre <- attr(Y.r, "scaled:center")
            CCA$Xbar <- Xbar
            Xbar <- qr.resid(Q, Xbar)
        } else {
            CCA <- list(eig = 0, rank = rank, qrank = qrank, tot.chi = 0,
                        QR = Q, Xbar = Xbar)
            u <- matrix(0, nrow=nrow(sol$u), ncol=0)
            v <- matrix(0, nrow=nrow(sol$v), ncol=0)
            CCA$u <- CCA$wa <- u
            CCA$v <- v
            CCA$biplot <- matrix(0, 0, 0)
            CCA$alias <- colnames(Y.r)
        }
    }
    Q <- qr(Xbar)
    sol <- svd(Xbar)
    sol$d <- sol$d/sqrt(NR)
    ax.names <- paste("PC", 1:length(sol$d), sep = "")
    colnames(sol$u) <- ax.names
    colnames(sol$v) <- ax.names
    names(sol$d) <- ax.names
    rownames(sol$u) <- rownames(X)
    rownames(sol$v) <- colnames(X)
    rank <- min(Q$rank, sum(sol$d > (sol$d[1L] * ZERO)))
    if (rank) {
        CA <- list(eig = (sol$d[1:rank]^2))
        CA$u <- as.matrix(sol$u)[, 1:rank, drop = FALSE]
        CA$v <- as.matrix(sol$v)[, 1:rank, drop = FALSE]
        CA$rank <- rank
        CA$tot.chi <- sum(CA$eig)
        CA$Xbar <- Xbar
    } else {   # zero rank: no residual component
        CA <- list(eig = 0, rank = rank, tot.chi = 0,
                   Xbar = Xbar)
        CA$u <- matrix(0, nrow(sol$u), 0)
        CA$v <- matrix(0, nrow(sol$v), 0)
    }
    call <- match.call()
    call[[1]] <- as.name("rda")
    sol <- list(call = call, grand.total = NA, rowsum = NA, colsum = SD, 
                tot.chi = tot.chi, pCCA = pCCA, CCA = CCA, CA = CA)
    sol$method <- "rda"
    sol$inertia <- if (scale) 
        "correlations"
    else "variance"
    class(sol) <- c("rda", "cca")
    sol
}
