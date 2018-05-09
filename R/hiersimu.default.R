hiersimu.default <-
function(y, x, FUN, location = c("mean", "median"),
         relative = FALSE, drop.highest = FALSE, nsimul=99,
         method = "r2dtable", ...)
{
    ## evaluate formula
    lhs <- as.matrix(y)
    if (missing(x))
        x <- cbind(level_1=seq_len(nrow(lhs)),
            leve_2=rep(1, nrow(lhs)))
    rhs <- data.frame(x)
    rhs[] <- lapply(rhs, as.factor)
    rhs[] <- lapply(rhs, droplevels, exclude = NA)
    nlevs <- ncol(rhs)
    if (is.null(colnames(rhs)))
        colnames(rhs) <- paste("level", 1:nlevs, sep="_")
    tlab <- colnames(rhs)

    ## check proper design of the model frame
    l1 <- sapply(rhs, function(z) length(unique(z)))
    if (!any(sapply(2:nlevs, function(z) l1[z] <= l1[z-1])))
        stop("number of levels are inappropriate, check sequence")
    rval <- list()
    rval[[1]] <- rhs[,nlevs]
    nCol <- nlevs - 1
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
        stop("levels are not perfectly nested")

    ## aggregate response matrix
    fullgamma <-if (nlevels(rhs[,nlevs]) == 1)
        TRUE else FALSE
    if (fullgamma && drop.highest)
        nlevs <- nlevs - 1
    if (nlevs == 1 && relative)
        stop("'relative=FALSE' makes no sense with one level")
    ftmp <- vector("list", nlevs)
    for (i in 1:nlevs) {
        ftmp[[i]] <- as.formula(paste("~", tlab[i], "- 1"))
    }

    ## is there burnin/thin in ... ?
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
    call <- match.call()
    call[[1]] <- as.name("hiersimu")
    attr(sim, "call") <- call
    attr(sim, "FUN") <- FUN
    attr(sim, "location") <- location
    attr(sim, "n.levels") <- nlevs
    attr(sim, "terms") <- tlab
    attr(sim, "model") <- rhs
    class(sim) <- c("hiersimu", class(sim))
    sim
}
