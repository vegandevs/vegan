adipart <-
function(formula, data, index=c("richness", "shannon", "simpson"),
    weights=c("unif", "prop"), relative = FALSE, nsimul=99, ...)
{
    ## evaluate formula
    lhs <- formula[[2]]
    if (missing(data))
        data <- parent.frame()
    lhs <- as.matrix(eval(lhs, data))
    formula[[2]] <- NULL
    rhs <- model.frame(formula, data, drop.unused.levels = TRUE)
    tlab <- attr(attr(rhs, "terms"), "term.labels")
    nlevs <- length(tlab)
    if (nlevs < 2)
        stop("provide at least two level hierarchy")
    if (any(rowSums(lhs) == 0))
        stop("data matrix contains empty rows")
    if (any(lhs < 0))
        stop("data matrix contains negative entries")

    ## part check proper design of the model frame
    noint <- attr(attr(attr(rhs, "terms"), "factors"), "dimnames")[[1]]
    int <- attr(attr(attr(rhs, "terms"), "factors"), "dimnames")[[2]]
    if (!identical(noint, int))
        stop("interactions are not allowed in formula")
    if (!all(attr(attr(rhs, "terms"), "dataClasses") == "factor"))
        stop("all right hand side variables in formula must be factors")
    l1 <- sapply(rhs, function(z) length(unique(z)))
    if (!any(sapply(2:nlevs, function(z) l1[z] <= l1[z-1])))
        stop("number of levels are inapropriate, check sequence")
    rval <- list()
    rval[[1]] <- as.factor(rhs[,nlevs])
    rval[[1]] <- rval[[1]][drop = TRUE]
    nCol <- nlevs - 1
    for (i in 2:nlevs) {
        rval[[i]] <- interaction(rhs[,nCol], rval[[(i-1)]], drop=TRUE)
        nCol <- nCol - 1
    }
    rval <- as.data.frame(rval[rev(1:length(rval))])
    l2 <- sapply(rval, function(z) length(unique(z)))
    if (any(l1 != l2))
        warning("levels are not perfectly nested")

    ## aggregate response matrix
    fullgamma <-if (nlevels(rhs[,nlevs]) == 1)
        TRUE else FALSE
    ftmp <- vector("list", nlevs)
    for (i in 1:nlevs) {
        ftmp[[i]] <- as.formula(paste("~", tlab[i], "- 1"))
    }

    ## is there a method/burnin/thin in ... ?
    method <- if (is.null(list(...)$method))
        "r2dtable" else list(...)$method
    burnin <- if (is.null(list(...)$burnin))
        0 else list(...)$burnin
    thin <- if (is.null(list(...)$thin))
        1 else list(...)$thin
    base <- if (is.null(list(...)$base))
        exp(1) else list(...)$base

    ## evaluate other arguments
    index <- match.arg(index)
    weights <- match.arg(weights)
    switch(index,
        "richness" = {
            divfun <- function(x) apply(x > 0, 1, sum)},
        "shannon" = {
            divfun <- function(x) diversity(x, index = "shannon", MARGIN = 1, base=base)},
        "simpson" = {
            divfun <- function(x) diversity(x, index = "simpson", MARGIN = 1)})
    sumMatr <- sum(lhs)

    ## this is the function passed to oecosimu
    wdivfun <- function(x) {
        if (fullgamma) {
            tmp <- lapply(1:(nlevs-1), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
            tmp[[nlevs]] <- matrix(colSums(lhs), nrow = 1, ncol = ncol(lhs))
        } else {
            tmp <- lapply(1:nlevs, function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
        }
        ## weights will change in oecosimu thus need to be recalculated
        if (weights == "prop")
            wt <- lapply(1:nlevs, function(i) apply(tmp[[i]], 1, function(z) sum(z) / sumMatr))
            else wt <- lapply(1:nlevs, function(i) rep(1 / NROW(tmp[[i]]), NROW(tmp[[i]])))
        a <- sapply(1:nlevs, function(i) sum(divfun(tmp[[i]]) * wt[[i]]))
        if (relative)
            a <- a / a[length(a)]
        b <- sapply(2:nlevs, function(i) a[i] - a[(i-1)])
        c(a, b)
    }
    if (nsimul > 0) {
        sim <- oecosimu(lhs, wdivfun, method = method, nsimul=nsimul,
            burnin=burnin, thin=thin)
        } else {
            sim <- wdivfun(lhs)
            tmp <- rep(NA, length(sim))
            sim <- list(statistic = sim,
                oecosimu = list(z = tmp, pval = tmp, method = NA, statistic = sim))
        }
    nam <- c(paste("alpha", 1:(nlevs-1), sep="."), "gamma",
        paste("beta", 1:(nlevs-1), sep="."))
    names(sim$statistic) <- attr(sim$oecosimu$statistic, "names") <- nam
    attr(sim, "call") <- match.call()
    attr(sim, "index") <- index
    attr(sim, "weights") <- weights
    attr(sim, "n.levels") <- nlevs
    attr(sim, "terms") <- tlab
    attr(sim, "model") <- rhs
    class(sim) <- c("adipart", "list")
    sim
}

