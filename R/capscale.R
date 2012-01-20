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
    X <- eval(formula[[2]], envir=parent.frame(),
              enclos = environment(formula))
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
    fla <- update(formula, X ~ .)
    environment(fla) <- environment()
    d <- ordiParseFormula(fla,
                          if(is.data.frame(data) && !is.null(comm)) cbind(data, comm)
                          else data,
                          envdepth = 1, na.action = na.action,
                          subset = substitute(subset))
    ## ordiParseFormula subsets rows of dissimilarities: do the same
    ## for columns ('comm' is handled later)
    if (!is.null(d$subset))
        d$X <- d$X[, d$subset, drop = FALSE]
    ## Delete columns if rows were deleted due to missing values
    if (!is.null(d$na.action)) {
        d$X <- d$X[, -d$na.action, drop = FALSE]
    }
    X <- as.dist(d$X)
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
        X <- cmdscale(X, k = k, eig = TRUE, add = add)
        ## All eigenvalues *should* be positive, but see that they are
        if (getRversion() < "2.13.0")
            X$points <- X$points[, X$eig[-(k+1)] > 0]
        X$eig <- X$eig[X$eig > 0]
    }
    else
        X <- wcmdscale(X, eig = TRUE)
    if (is.null(rownames(X$points))) 
        rownames(X$points) <- nm
    X$points <- adjust * X$points
    if (adjust == 1)
        X$eig <- X$eig/k
    sol <- rda.default(X$points, d$Y, d$Z, ...)
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
        sol$tot.chi <- sol$tot.chi + sol$CA$imaginary.chi
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
            sol$CCA$v.eig <- t(comm) %*% sol$CCA$u/sqrt(k)
            sol$CCA$v <- sweep(sol$CCA$v.eig, 2, sqrt(sol$CCA$eig), 
                               "/")
            comm <- qr.resid(sol$CCA$QR, comm)
        }
        if (!is.null(sol$CA) && sol$CA$rank > 0) {
            sol$CA$v.eig <- t(comm) %*% sol$CA$u/sqrt(k)
            sol$CA$v <- sweep(sol$CA$v.eig, 2, sqrt(sol$CA$eig[1:poseig]), 
                              "/")
        }
    } else {
        ## input data were dissimilarities, and no 'comm' defined:
        ## species scores make no sense and are made NA
        sol$CA$v.eig[] <- sol$CA$v[] <- NA
        if (!is.null(sol$CCA))
            sol$CCA$v.eig[] <- sol$CCA$v[] <- NA
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
    if (add)
        sol$ac <- X$ac
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
