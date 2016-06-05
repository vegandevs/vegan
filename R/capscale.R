`capscale` <-
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
    ## see if user supplied dissimilarities as a matrix
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
    inertia <- attr(X, "method")
    if (is.null(inertia))
        inertia <- "unknown"
    inertia <- paste(toupper(substr(inertia, 1, 1)),
                     substring(inertia,  2), sep = "")
    inertia <- paste(inertia, "distance")
    if (!sqrt.dist)
        inertia <- paste("squared", inertia)
    ## postpone info on euclidification till we have done so

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
    ## wcmdscale, optionally with additive adjustment
    X <- wcmdscale(X, x.ret = TRUE, add = add)
    ## this may have been euclidified: update inertia
    if (!is.na(X$ac) && X$ac > sqrt(.Machine$double.eps))
        inertia <- paste(paste0(toupper(substring(X$add, 1, 1)),
                                substring(X$add, 2)),
                         "adjusted", inertia)
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
    sol <- rda.default(X$points, d$Y, d$Z, ...)
    ## Get components of inertia with negative eigenvalues following
    ## McArdle & Anderson (2001), section "Theory". G is their
    ## double-centred Gower matrix, but instead of hat matrix, we use
    ## use QR decomposition to get the components of inertia.
    hasNegEig <- any(X$eig < 0)
    G <- -X$x/2
    if (adjust == 1)
        G <- G/k
    if (hasNegEig)
        sol$real.tot.chi <- sol$tot.chi
    sol$tot.chi <- sum(diag(G))
    if (!is.null(sol$pCCA)) {
        sol$pCCA$G <- G
        if (hasNegEig) {
            sol$pCCA$real.tot.chi <- sol$pCCA$tot.chi
            sol$pCCA$tot.chi <- sum(diag(qr.fitted(sol$pCCA$QR, G)))
        }
        G <- qr.resid(sol$pCCA$QR, t(qr.resid(sol$pCCA$QR, G)))
    }
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
        sol$CCA$G <- G
        if (hasNegEig) {
            sol$CCA$real.tot.chi <- sol$CCA$tot.chi
            sol$CCA$tot.chi <- sum(diag(qr.fitted(sol$CCA$QR, G)))
        }
    }
    if (hasNegEig) {
        sol$CA$real.tot.chi <- sol$CA$tot.chi
        if (!is.null(sol$CA) && !is.null(sol$CCA$QR))
            sol$CA$tot.chi <- sum(diag(qr.resid(sol$CCA$QR, G)))
        else
            sol$CA$tot.chi <- sum(diag(G))
    }
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
        colnames(sol$CCA$u) <- colnames(sol$CCA$biplot) <- names(sol$CCA$eig) <-
            colnames(sol$CCA$wa) <- colnames(sol$CCA$v) <-
                paste("CAP", 1:ncol(sol$CCA$u), sep = "")
    }
    if (!is.null(sol$CA) && sol$CA$rank > 0) {
        colnames(sol$CA$u) <- names(sol$CA$eig) <- colnames(sol$CA$v) <-
            paste("MDS", 1:ncol(sol$CA$u), sep = "")
    }
    ## update for negative eigenvalues
    poseig <- length(sol$CA$eig)
    if (any(X$eig < 0)) {
        negax <- X$eig[X$eig < 0]
        sol$CA$imaginary.chi <- sum(negax)
        sol$CA$imaginary.rank <- length(negax)
        sol$CA$imaginary.u.eig <- X$negaxes
    }
    if (!is.null(comm)) {
        comm <- scale(comm, center = TRUE, scale = FALSE)
        sol$colsum <- apply(comm, 2, sd)
        ## take a 'subset' of the community after scale()
        if (!is.null(d$subset))
            comm <- comm[d$subset, , drop = FALSE]
        ## NA action after 'subset'
        if (!is.null(d$na.action))
            comm <- comm[-d$na.action, , drop = FALSE]
        if (!is.null(sol$pCCA) && sol$pCCA$rank > 0) 
            comm <- qr.resid(sol$pCCA$QR, comm)
        if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
            v.eig <- t(comm) %*% sol$CCA$u/sqrt(k)
            sol$CCA$v <- decostand(v.eig, "normalize", MARGIN = 2)
            comm <- qr.resid(sol$CCA$QR, comm)
        }
        if (!is.null(sol$CA) && sol$CA$rank > 0) {
            v.eig <- t(comm) %*% sol$CA$u/sqrt(k)
            sol$CA$v <- decostand(v.eig, "normalize", MARGIN = 2)
        }
    } else {
        ## input data were dissimilarities, and no 'comm' defined:
        ## species scores make no sense and are made NA
        sol$CA$v[] <- NA
        if (!is.null(sol$CCA))
            sol$CCA$v[] <- NA
        sol$colsum <- NA
    }
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) 
        sol$CCA$centroids <- centroids.cca(sol$CCA$wa, d$modelframe)
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
    sol$method <- "capscale"
    sol$sqrt.dist <- sqrt.dist
    if (!is.na(X$ac) && X$ac > 0) {
        sol$ac <- X$ac
        sol$add <- X$add
    }
    sol$adjust <- adjust
    sol$inertia <- inertia
    if (metaMDSdist)
        sol$metaMDSdist <- commname
    sol$subset <- d$subset
    sol$na.action <- d$na.action
    class(sol) <- c("capscale", class(sol))
    if (!is.null(sol$na.action))
        sol <- ordiNAexclude(sol, d$excluded)
    sol
}
