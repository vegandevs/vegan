multipart <-
function(formula, data, index=c("renyi", "tsallis"), scales = 1,
    nsimul=99, control, ...)
{
    if (length(scales) > 1)
        stop("length of 'scales' must be 1")
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
    if (!fullgamma)
        warning("gamma diversity value might be meaningless (untested feature)")
    ftmp <- vector("list", nlevs)
    for (i in 1:nlevs) {
        ftmp[[i]] <- as.formula(paste("~", tlab[i], "- 1"))
    }

    ## evaluate other arguments
    index <- match.arg(index)
    if (missing(control))
        control <- permat.control()
    divfun <- switch(index,
        "renyi" = function(x) renyi(x, scales=scales, hill = TRUE),
        "tsallis" = function(x) tsallis(x, scales=scales, hill = TRUE))

    ## this is the function passed to oecosimu
    wdivfun <- function(x) {
        if (fullgamma) {
            tmp <- lapply(1:(nlevs-1), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
            tmp[[nlevs]] <- matrix(colSums(lhs), nrow = 1, ncol = ncol(lhs))
        } else {
            tmp <- lapply(1:nlevs, function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
        }
        a <- sapply(1:nlevs, function(i) mean(divfun(tmp[[i]]), na.rm=TRUE))
        G <- a[nlevs]
        b <- sapply(1:(nlevs-1), function(i) G / a[i])
        c(a, b)
    }
    sim <- oecosimu(lhs, wdivfun, method = "permat", nsimul=nsimul,
        burnin=control$burnin, thin=control$thin, control=control)
    nam <- c(paste("alpha", 1:(nlevs-1), sep="."), "gamma",
        paste("beta", 1:(nlevs-1), sep="."))
    names(sim$statistic) <- attr(sim$oecosimu$statistic, "names") <- nam
    attr(sim, "call") <- match.call()
    attr(sim, "index") <- index
    attr(sim, "scales") <- scales
    attr(sim, "global") <- TRUE
    attr(sim, "n.levels") <- nlevs
    attr(sim, "terms") <- tlab
    attr(sim, "model") <- rhs
    class(sim) <- c("multipart", "list")
    sim
}
