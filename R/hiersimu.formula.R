hiersimu.formula <-
function(formula, data, FUN, location = c("mean", "median"),
relative = FALSE, drop.highest = FALSE, nsimul=99, ...)
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

    ## part check proper design of the model frame
    noint <- attr(attr(attr(rhs, "terms"), "factors"), "dimnames")[[1]]
    int <- attr(attr(attr(rhs, "terms"), "factors"), "dimnames")[[2]]
    if (!identical(noint, int))
        stop("interactions are not allowed in formula")
    if (!all(attr(attr(rhs, "terms"), "dataClasses") == "factor"))
        stop("all right hand side variables in formula must be factors")
    l1 <- sapply(rhs, function(z) length(unique(z)))
    if (nlevs > 1 && !any(sapply(2:nlevs, function(z) l1[z] <= l1[z-1])))
        stop("number of levels are inapropriate, check sequence")
    rval <- list()
    rval[[1]] <- as.factor(rhs[,nlevs])
    rval[[1]] <- rval[[1]][drop = TRUE]
    if (nlevs > 1) {
        nCol <- nlevs - 1
        for (i in 2:nlevs) {
            rval[[i]] <- interaction(rhs[,nCol], rval[[(i-1)]], drop=TRUE)
            nCol <- nCol - 1
        }
    }
    rval <- as.data.frame(rval[rev(1:length(rval))])
    l2 <- sapply(rval, function(z) length(unique(z)))
    if (any(l1 != l2))
        warning("levels are not perfectly nested")

    ## aggregate response matrix
    fullgamma <-if (nlevels(rhs[,nlevs]) == 1)
        TRUE else FALSE
    if (fullgamma && drop.highest)
        nlevs <- nlevs - 1
    if (nlevs == 1 && relative)
        stop("'relative=FALSE' makes no sense with 1 level")
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

    ## evaluate other arguments
    if (!is.function(FUN))
        stop("'FUN' must be a function")
    location <- match.arg(location)
    aggrFUN <- switch(location,
        "mean" = mean,
        "median" = median)

    ## this is the function passed to oecosimu
    evalFUN <- function(x) {
        if (fullgamma && !drop.highest) {
            tmp <- lapply(1:(nlevs-1), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
            tmp[[nlevs]] <- matrix(colSums(x), nrow = 1, ncol = ncol(x))
        } else {
            tmp <- lapply(1:nlevs, function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
        }
        a <- sapply(1:nlevs, function(i) aggrFUN(FUN(tmp[[i]]))) # dots removed from FUN
        if (relative)
            a <- a / a[length(a)]
        a
    }

    ## processing oecosimu results
    sim <- oecosimu(lhs, evalFUN, method = method, nsimul=nsimul,
        burnin=burnin, thin=thin)
#    nam <- paste("level", 1:nlevs, sep=".")
#    names(sim$statistic) <- attr(sim$oecosimu$statistic, "names") <- nam
    names(sim$statistic) <- attr(sim$oecosimu$statistic, "names") <- tlab[1:nlevs]
    attr(sim, "call") <- match.call()
    attr(sim, "FUN") <- FUN
    attr(sim, "location") <- location
    attr(sim, "n.levels") <- nlevs
    attr(sim, "terms") <- tlab
    attr(sim, "model") <- rhs
    class(sim) <- c("hiersimu", "list")
    sim
}
