`dbrda` <-
    function (formula, data, distance = "euclidean",
              sqrt.dist = FALSE,  add = FALSE, dfun = vegdist,
              metaMDSdist = FALSE, na.action = na.fail,
              subset = NULL, ...) 
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
    if ((is.matrix(X) || is.data.frame(X)) &&
               isSymmetric(unname(as.matrix(X))))
        X <- as.dist(X)
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
    ## get the name of the inertia
    inertia <- attr(X, "method")
    if (is.null(inertia))
        inertia <- "unknown"
    inertia <- paste(toupper(substr(inertia, 1, 1)),
                     substring(inertia, 2), sep = "")
    inertia <- paste(inertia, "distance")

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
    X <- as.matrix(X)
    k <- NROW(X) - 1
    ## sqrt & add adjustments
    if (sqrt.dist)
        X <- sqrt(X)
    if (is.logical(add) && isTRUE(add))
        add <- "lingoes"
    if (is.character(add)) {
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- addLingoes(X)
            X <- sqrt(X^2 + 2 * ac)
        } else if (add == "cailliez") {
            ac <- addCailliez(X)
            X <- X + ac
        }
        diag(X) <- 0
    } else {
        ac <- 0
    }
    ## update the name of the inertia
    if (!sqrt.dist)
        inertia <- paste("squared", inertia)
    if (ac > sqrt(.Machine$double.eps))
        inertia <- paste(paste0(toupper(substring(add, 1, 1)),
                              substring(add, 2)), "adjusted", inertia)
    if (max(X) >= 4 + .Machine$double.eps) {
        inertia <- paste("mean", inertia)
        adjust <- 1
    }
    else {
        adjust <- sqrt(k)
    }
    nm <- attr(X, "Labels")    
    ## Get components of inertia with negative eigenvalues following
    ## McArdle & Anderson (2001), section "Theory". G is their
    ## double-centred Gower matrix, but instead of hat matrix, we use
    ## QR decomposition to get the components of inertia.
    G <- -GowerDblcen(X^2)/2
    if (adjust == 1)
        G <- G/k
    ## Solution: this shows the algorithmic steps
    tot.chi <- sum(diag(G))
    pCCA <- CCA <-  CA <- NULL
    ## pCCA
    if (!is.null(d$Z)) {
        d$Z <- scale(d$Z, scale = FALSE)
        Q <- qr(d$Z, tol = 1e-6)
        HGH <- qr.fitted(Q, t(qr.fitted(Q, G)))
        pCCA <- list(rank = Q$rank, tot.chi = sum(diag(HGH)),
                     QR = Q, Fit = HGH,
                     envcentre = attr(d$Z, "scaled:center"),
                     G = G)
        G <- qr.resid(Q, t(qr.resid(Q, G)))
    }
    ## CCA
    if (!is.null(d$Y)) {
        d$Y <- scale(d$Y, scale = FALSE) 
        Q <- qr(cbind(d$Z, d$Y), tol = 1e-6)
        HGH <- qr.fitted(Q, t(qr.fitted(Q, G)))
        e <- eigen(HGH, symmetric = TRUE)
        nz <- abs(e$values) > EPS
        if (any(nz)) {
            e$values <- e$values[nz]
            e$vectors <- e$vectors[, nz, drop = FALSE]
            pos <- e$values > 0
            if (any(e$values < 0)) {
                imaginary.u <- e$vectors[, !pos, drop = FALSE]
                e$vectors <- e$vectors[, pos, drop = FALSE]
            } else {
                imaginary.u <- NULL
            }
            wa <- G %*% e$vectors %*% diag(1/e$values[pos], sum(pos))
            v <- matrix(NA, ncol = ncol(wa))
            oo <- Q$pivot[seq_len(Q$rank)]
            rank <- Q$rank
            if (!is.null(pCCA)) {
                oo <- oo[-seq_len(pCCA$rank)] - ncol(d$Z)
                rank <- rank - pCCA$rank
            }
            CCA <- list(eig = e$values,
                        u = e$vectors,
                        imaginary.u = imaginary.u,
                        poseig = sum(pos),
                        v = v, wa = wa,
                        alias =  if (rank < ncol(d$Y))
                                     colnames(d$Y)[-oo],
                        biplot = cor(d$Y[,oo, drop=FALSE], e$vectors),
                        qrank = rank, rank = rank,
                        tot.chi = sum(diag(HGH)),
                        QR = Q,
                        envcentre = attr(d$Y, "scaled:center"),
                        Xbar = NA, G = G)
        } else {
            CCA <- NULL
        }
        G <- qr.resid(Q, t(qr.resid(Q, G)))
    }
    ## CA
    e <- eigen(G, symmetric = TRUE)
    nz <- abs(e$values) > EPS # positively or negatively non-zero
    if (any(nz)) {
        e$values <- e$values[nz]
        e$vectors <- e$vectors[, nz, drop = FALSE]
        if (any(e$values < 0)) {
            imaginary.u <- e$vectors[, e$values < 0, drop = FALSE]
            e$vectors <- e$vectors[, e$values > 0, drop = FALSE]
        } else {
            imaginary.u <- NULL
        }
        v <- matrix(NA, ncol = ncol(e$vectors))
        CA <- list(eig = e$values,
                   u = e$vectors,
                   imaginary.u = imaginary.u,
                   poseig = sum(e$values > 0),
                   v = v,
                   rank = sum(nz),
                   tot.chi = sum(diag(G)),
                   Xbar = NA, G = G)
    } else {
        CA <- NULL
    }
    ## output
    sol <- list(tot.chi = tot.chi, pCCA = pCCA, CCA = CCA, CA = CA)
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
        colnames(sol$CCA$u) <-
            colnames(sol$CCA$wa) <-
            colnames(sol$CCA$v) <-
            names(sol$CCA$eig) <-
                paste("dbRDA", seq_len(ncol(sol$CCA$u)), sep = "")
        colnames(sol$CCA$biplot) <-
            names(sol$CCA$eig)[sol$CCA$eig > 0]
        rownames(sol$CCA$u) <- rownames(d$X)
        if (!is.null(sol$CCA$imaginary.u)) {
            negax <- sol$CCA$eig < 0
            negnm <- paste0("idbRDA", seq_len(sum(negax)))
            names(sol$CCA$eig)[negax] <- negnm
            colnames(sol$CCA$imaginary.u) <- negnm
            rownames(sol$CCA$imaginary.u) <- rownames(d$X)
        }
    }
    if (!is.null(sol$CA) && sol$CA$rank > 0) {
        colnames(sol$CA$u) <- colnames(sol$CA$v) <- names(sol$CA$eig) <-
            paste("MDS", seq_len(ncol(sol$CA$u)), sep = "")
        rownames(sol$CA$u) <- rownames(d$X)
        if (!is.null(sol$CA$imaginary.u)) {
            negax <- sol$CA$eig < 0
            negnm <- paste0("iMDS", seq_len(sum(negax)))
            names(sol$CA$eig)[negax] <- negnm
            colnames(sol$CA$imaginary.u) <- negnm
            rownames(sol$CA$imaginary.u) <- rownames(d$X)
        }
    }

    sol$colsum <- NA
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) 
        sol$CCA$centroids <-
            centroids.cca(sol$CCA$u, d$modelframe)
    if (!is.null(sol$CCA$alias)) 
        sol$CCA$centroids <- unique(sol$CCA$centroids)
    if (!is.null(sol$CCA$centroids)) {
        rs <- rowSums(sol$CCA$centroids^2)
        sol$CCA$centroids <- sol$CCA$centroids[rs > 1e-04, , 
                                               drop = FALSE]
        if (nrow(sol$CCA$centroids) == 0)
            sol$CCA$centroids <- NULL
    }
    sol$call <- match.call()
    sol$terms <- terms(formula, "Condition", data = data)
    sol$terminfo <- ordiTerminfo(d, data)
    sol$call$formula <- formula(d$terms, width.cutoff = 500)
    sol$call$formula[[2]] <- formula[[2]]
    sol$method <- "dbrda"
    sol$sqrt.dist <- sqrt.dist
    if (!is.na(ac) && ac > 0) {
        sol$ac <- ac
        sol$add <- add
    }
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
