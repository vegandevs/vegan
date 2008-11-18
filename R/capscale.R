`capscale` <-
    function (formula, data, distance = "euclidean", comm = NULL, 
              add = FALSE, dfun = vegdist, metaMDSdist = FALSE, ...) 
{
    if (!inherits(formula, "formula")) 
        stop("Needs a model formula")
    if (missing(data)) {
        data <- parent.frame()
    }
    else {
        data <- ordiGetData(match.call(), environment(formula))
    }
    formula <- formula(terms(formula, data = data))
    X <- formula[[2]]
    X <- eval(X, environment(formula))
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
    inertia <- paste("squared", inertia, "distance")
    if (add) 
        inertia <- paste(inertia, "(euclidified)")
    k <- attr(X, "Size") - 1
    if (max(X) >= 4 + .Machine$double.eps) {
        inertia <- paste("mean", inertia)
        adjust <- 1
    }
    else {
        adjust <- k
    }
    nm <- attr(X, "Labels")
    X <- cmdscale(X, k = k, eig = TRUE, add = add)
    if (is.null(rownames(X$points))) 
        rownames(X$points) <- nm
    X$points <- adjust * X$points
    neig <- min(which(X$eig < 0) - 1, k)
    sol <- X$points[, 1:neig]
    fla <- update(formula, sol ~ .)
    environment(fla) <- environment()
    d <- ordiParseFormula(fla, data, envdepth = 1)
    sol <- rda.default(d$X, d$Y, d$Z, ...)
    sol$tot.chi <- sol$tot.chi
    if (!is.null(sol$CCA)) {
        colnames(sol$CCA$u) <- colnames(sol$CCA$biplot) <- names(sol$CCA$eig) <-
            colnames(sol$CCA$wa) <- colnames(sol$CCA$v) <-
                paste("CAP", 1:ncol(sol$CCA$u), sep = "")
    }
    if (!is.null(sol$CA)) {
        colnames(sol$CA$u) <- names(sol$CA$eig) <- colnames(sol$CA$v) <-
            paste("MDS", 1:ncol(sol$CA$u), sep = "")
    }
    if (!is.null(comm)) {
        comm <- scale(comm, center = TRUE, scale = FALSE)
        sol$colsum <- sd(comm)
        if (!is.null(sol$pCCA)) 
            comm <- qr.resid(sol$pCCA$QR, comm)
        if (!is.null(sol$CCA)) {
            sol$CCA$v.eig <- t(comm) %*% sol$CCA$u/sqrt(k)
            sol$CCA$v <- sweep(sol$CCA$v.eig, 2, sqrt(sol$CCA$eig), 
                               "/")
            comm <- qr.resid(sol$CCA$QR, comm)
        }
        if (!is.null(sol$CA)) {
            sol$CA$v.eig <- t(comm) %*% sol$CA$u/sqrt(k)
            sol$CA$v <- sweep(sol$CA$v.eig, 2, sqrt(sol$CA$eig), 
                              "/")
        }
    }
    if (!is.null(sol$CCA)) 
        sol$CCA$centroids <- centroids.cca(sol$CCA$wa, d$modelframe)
    if (!is.null(sol$CCA$alias)) 
        sol$CCA$centroids <- unique(sol$CCA$centroids)
    if (!is.null(sol$CCA$centroids)) {
        rs <- rowSums(sol$CCA$centroids^2)
        sol$CCA$centroids <- sol$CCA$centroids[rs > 1e-04, , 
                                               drop = FALSE]
    }
    sol$call <- match.call()
    sol$terms <- terms(formula, "Condition", data = data)
    sol$terminfo <- ordiTerminfo(d, data)
    sol$call$formula <- formula(d$terms, width.cutoff = 500)
    sol$call$formula[[2]] <- formula[[2]]
    sol$method <- "capscale"
    sol$inertia <- inertia
    if (metaMDSdist)
        sol$metaMDSdist <- commname
    class(sol) <- c("capscale", class(sol))
    sol
}
