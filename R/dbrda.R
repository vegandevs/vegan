`dbrda` <-
    function (formula, data, distance = "euclidean", sqrt.dist = FALSE,
              comm = NULL, add = FALSE, dfun = vegdist,
              metaMDSdist = FALSE, na.action = na.fail, subset = NULL, ...) 
{
    EPS <- sqrt(.Machine$double.eps)
    if (!inherits(formula, "formula")) 
        stop("Needs a model formula")
    if (missing(data)) {
        data <- parent.frame()
    }
    else {
        data <- ordiGetData(match.call(), environment(formula))
    }
    formula <- formula(terms(formula, data = data))
    ## The following line was eval'ed in environment(formula), but
    ## that made update() fail. Rethink the line if capscale() fails
    ## mysteriously at this point.
    X <- eval(formula[[2]], envir=environment(formula),
              enclos = globalenv())
    if (!inherits(X, "dist")) {
        comm <- X
        dfun <- match.fun(dfun)
        if (metaMDSdist) {
            commname <- as.character(formula[[2]])
            X <- metaMDSdist(comm, distance = distance, zerodist = "ignore",
                             commname = commname, distfun = dfun, ...)
            commname <- attr(X, "commname")
            comm <- eval.parent(parse(text=commname))
        } else {
            X <- dfun(X, distance)
        }
    }
    inertia <- attr(X, "method")
    if (is.null(inertia))
        inertia <- "unknown"
    inertia <- paste(toupper(substr(inertia, 1, 1)), substr(inertia, 
                                                            2, 256), sep = "")
    inertia <- paste(inertia, "distance")
    if (!sqrt.dist)
        inertia <- paste("squared", inertia)
    if (add) 
        inertia <- paste(inertia, "(euclidified)")

    ## evaluate formula: ordiParseFormula will return dissimilarities
    ## as a symmetric square matrix (except that some rows may be
    ## deleted due to missing values)
    d <- ordiParseFormula(formula,
                          data,
                          na.action = na.action,
                          subset = substitute(subset))
    ## ordiParseFormula subsets rows of dissimilarities: do the same
    ## for columns ('comm' is handled later). ordiParseFormula
    ## returned the original data, but we use instead the potentially
    ## changed X and discard d$X.
    if (!is.null(d$subset)) {
        X <- as.matrix(X)[d$subset, d$subset, drop = FALSE]
    }
    ## Delete columns if rows were deleted due to missing values
    if (!is.null(d$na.action)) {
        X <- as.matrix(X)[-d$na.action, -d$na.action, drop = FALSE]
    }
    X <- as.dist(X)
    k <- attr(X, "Size") - 1
    if (sqrt.dist)
        X <- sqrt(X)
    if (max(X) >= 4 + .Machine$double.eps) {
        inertia <- paste("mean", inertia)
        adjust <- 1
    }
    else {
        adjust <- sqrt(k)
    }
    nm <- attr(X, "Labels")    
    ## cmdscale is only used if 'add = TRUE': it cannot properly
    ## handle negative eigenvalues and therefore we normally use
    ## wcmdscale. If we have 'add = TRUE' there will be no negative
    ## eigenvalues and this is not a problem.
    if (add) {
        X <- cmdscale(X, k = k, eig = TRUE, add = add, x.ret = TRUE)
        ## All eigenvalues *should* be positive, but see that they are
        X$eig <- X$eig[X$eig > 0]
    }
    else
        X <- wcmdscale(X, x.ret = TRUE)
    if (is.null(rownames(X$points))) 
        rownames(X$points) <- nm
    X$points <- adjust * X$points
    ## We adjust eigenvalues to variances, and simultaneously the
    ## possible negative axes must be adjusted similarly
    if (adjust == 1) {
        X$eig <- X$eig/k
        if (!is.null(X$negaxes))
            X$negaxes <- X$negaxes/sqrt(k)
    }
    ## Get components of inertia with negative eigenvalues following
    ## McArdle & Anderson (2001), section "Theory". G is their
    ## double-centred Gower matrix, but instead of hat matrix, we use
    ## QR decomposition to get the components of inertia.
    hasNegEig <- any(X$eig < 0)
    G <- -X$x/2
    if (adjust == 1)
        G <- G/k
    ## Solution: this shows the algorithmic steps
    tot.chi <- sum(diag(G))
    pCCA <- CCA <-  CA <- NULL
    ## pCCA
    if (!is.null(d$Z)) {
        d$Z <- scale(d$Z, scale = FALSE)
        Q <- qr(d$Z, tol = 1e-6)
        H <- tcrossprod(qr.Q(Q)[, seq_len(Q$rank), drop=FALSE])
        HGH <- H %*% G %*% H
        pCCA <- list(rank = Q$rank, tot.chi = sum(diag(HGH)), QR = Q,
                     Fit = HGH, envcentre = attr(d$Z, "scaled:center"),
                     G = G)
        G <- G - HGH
    }
    ## CCA
    if (!is.null(d$Y)) {
        d$Y <- scale(d$Y, scale = FALSE) 
        Q <- qr(cbind(d$Z, d$Y), tol = 1e-6)
        H <- tcrossprod(qr.Q(Q)[, seq_len(Q$rank), drop=FALSE])
        HGH <- H %*% G %*% H
        e <- eigen(HGH)
        nz <- abs(e$values) > EPS
        e$values <- e$values[nz]
        e$vectors <- e$vectors[, nz, drop = FALSE]
        pos <- e$values > 0
        oo <- Q$pivot[seq_len(Q$rank)]
        rank <- Q$rank
        if (!is.null(pCCA)) {
            oo <- oo[oo > pCCA$rank] - ncol(d$Z)
            rank <- rank - pCCA$rank
        }
        CCA <- list(eig = e$values,
                    u = e$vectors,
                    v = NA, wa = NA,
                    biplot = cor(d$Y[,oo, drop=FALSE],
                    e$vectors[, pos, drop=FALSE]),
                    qrank = rank, rank = rank,
                    tot.chi = sum(diag(HGH)),
                    QR = Q,
                    envcentre = attr(d$Y, "scaled:center"),
                    Xbar = G, G = G)
        G <- G - HGH
    }
    ## CA
    e <- eigen(G) 
    nz <- abs(e$values) > EPS
    CA <- list(eig = e$values[nz],
               u = e$vectors[, nz, drop = FALSE],
               v = NA,
               rank = sum(nz),
               tot.chi = sum(diag(G)),
               Xbar = G)
    ## output
    sol <- list(tot.chi = tot.chi, pCCA = pCCA, CCA = CCA, CA = CA)
    ##if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
    ##   colnames(sol$CCA$u) <- colnames(sol$CCA$biplot) <- names(sol$CCA$eig) <-
    ##      colnames(sol$CCA$wa) <- colnames(sol$CCA$v) <-
    ##         paste("CAP", 1:ncol(sol$CCA$u), sep = "")
    ##}
    ##if (!is.null(sol$CA) && sol$CA$rank > 0) {
    ##    colnames(sol$CA$u) <- names(sol$CA$eig) <- colnames(sol$CA$v) <-
    ##        paste("MDS", 1:ncol(sol$CA$u), sep = "")
    ##}

    ## input data were dissimilarities, and no 'comm' defined:
    ## species scores make no sense and are made NA

    sol$colsum <- NA

    ## if (!is.null(sol$CCA) && sol$CCA$rank > 0) 
    ##    sol$CCA$centroids <- centroids.cca(sol$CCA$wa, d$modelframe)
    ## if (!is.null(sol$CCA$alias)) 
    ##    sol$CCA$centroids <- unique(sol$CCA$centroids)
    ## if (!is.null(sol$CCA$centroids)) {
    ##    rs <- rowSums(sol$CCA$centroids^2)
    ##    sol$CCA$centroids <- sol$CCA$centroids[rs > 1e-04, , 
    ##                                           drop = FALSE]
    ##    if (nrow(sol$CCA$centroids) == 0)
    ##        sol$CCA$centroids <- NULL
    ##}
    sol$call <- match.call()
    sol$terms <- terms(formula, "Condition", data = data)
    sol$terminfo <- ordiTerminfo(d, data)
    sol$call$formula <- formula(d$terms, width.cutoff = 500)
    sol$call$formula[[2]] <- formula[[2]]
    sol$method <- "capscale"
    if (add)
        sol$ac <- X$ac
    sol$adjust <- adjust
    sol$inertia <- inertia
    if (metaMDSdist)
        sol$metaMDSdist <- commname
    sol$subset <- d$subset
    sol$na.action <- d$na.action
    class(sol) <- c("dbrda", "rda", "cca")
    if (!is.null(sol$na.action))
        sol <- ordiNAexclude(sol, d$excluded)
    sol
}
