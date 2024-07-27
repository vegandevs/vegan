`dbrda` <-
    function (formula, data, distance = "euclidean",
              sqrt.dist = FALSE,  add = FALSE, dfun = vegdist,
              metaMDSdist = FALSE, na.action = na.fail,
              subset = NULL, ...)
{
    if (!inherits(formula, "formula"))
        stop("needs a model formula")
    if (missing(data)) {
        data <- parent.frame()
    }
    else {
        data <- eval(match.call()$data, environment(formula),
                     enclos = .GlobalEnv)
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
            X <- dfun(X, distance, ...)
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
                          subset = substitute(subset),
                          X = X)
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
        adjust <- sqrt(k)
        X <- X/adjust
    }
    else {
        adjust <- 1
    }
    ## Get components of inertia with negative eigenvalues following
    ## McArdle & Anderson (2001), section "Theory". G is their
    ## double-centred Gower matrix, but instead of hat matrix, we use
    ## QR decomposition to get the components of inertia.
    sol <- ordConstrained(X, d$Y, d$Z, method = "dbrda")

    sol$colsum <- NA
    ## separate eigenvectors associated with negative eigenvalues from
    ## u into imaginary.u
    if (!is.null(sol$CCA) && sol$CCA$rank > sol$CCA$poseig) {
        if (sol$CCA$poseig > 0)
            sol$CCA$imaginary.u <- sol$CCA$u[, -seq_len(sol$CCA$poseig),
                                             drop = FALSE]
        else
            sol$CCA$imaginary.u <- sol$CCA$u
        sol$CCA$u <- sol$CCA$u[, seq_len(sol$CCA$poseig), drop = FALSE]
    }
    if (!is.null(sol$CA) && sol$CA$rank > sol$CA$poseig) {
        sol$CA$imaginary.u <- sol$CA$u[, -seq_len(sol$CA$poseig),
                                       drop = FALSE]
        sol$CA$u <- sol$CA$u[, seq_len(sol$CA$poseig), drop = FALSE]
    }
    sol$CCA$centroids <- getCentroids(sol, d$modelframe)

    sol$call <- match.call()
    sol$terms <- terms(formula, "Condition", data = data)
    sol$terminfo <- ordiTerminfo(d, data)
    sol$call$formula <- formula(d$terms, width.cutoff = 500)
    sol$call$formula[[2]] <- formula[[2]]
    sol$sqrt.dist <- sqrt.dist
    if (!is.na(ac) && ac > 0) {
        sol$ac <- ac
        sol$add <- add
    }
    sol$adjust <- adjust
    sol$inertia <- inertia
    if (metaMDSdist)
        sol$metaMDSdist <- commname
    if (!is.null(d$subset))
        sol$subset <- d$subset
    if (!is.null(d$na.action)) {
        sol$na.action <- d$na.action
        ## dbrda cannot add WA scores in na.exclude, and the following
        ## does nothing except adds residuals.zombie
        sol <- ordiNAexclude(sol, d$excluded)
    }
    class(sol) <- c("dbrda", "rda", "cca")
    sol
}
