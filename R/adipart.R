adipart <-
function(matr, strata, index=c("richness", "shannon", "simpson"),
    weights=c("unif", "prop"), nsimul=99, control, ...)
{

    ## internal function, maybe later gets out
    nested.factor <-
    function(x) {
        x <- data.frame(x)
        nc <- NCOL(x)
        if (nc < 2)
            stop("number of columns must at least 2")
        nr <- NROW(x)
        l1 <- sapply(x, function(z) length(unique(z)))
        if (!any(sapply(2:nc, function(z) l1[z] <= l1[z-1])))
            stop("number of levels are inapropriate")
        rval <- list()
        rval[[1]] <- as.factor(x[,nc])
        rval[[1]] <- rval[[1]][drop = TRUE]
        ncol <- nc - 1
        for (i in 2:nc) {
            rval[[i]] <- interaction(x[,ncol], rval[[(i-1)]], drop=TRUE)
            ncol <- ncol - 1
        }
        rval <- as.data.frame(rval[rev(1:length(rval))])
        colnames(rval) <- paste("x", 1:nc, sep="")
        l2 <- sapply(rval, function(z) length(unique(z)))
        if (any(l1 != l2))
            warning("levels are not perfectly nested")
        rval
    }

    index <- match.arg(index)
    weights <- match.arg(weights)
    if (missing(control))
        control <- permat.control()
    switch(index,
        "richness" = {
            divfun <- function(x) apply(x > 0, 1, sum)},
        "shannon" = {
            divfun <- function(x) diversity(x, index = "shannon", MARGIN = 1, ...)},
        "simpson" = {
            divfun <- function(x) diversity(x, index = "simpson", MARGIN = 1)})
    strata <- nested.factor(strata)
    nl <- NCOL(strata)
#    seed <- trunc(runif(1, 1000, 9999))
    ## this is the function passed to oecosimu
    wdivfun <- function(x) {
        tmp <- lapply(1:nl, function(i) aggregate(x, list(strata[,i]), sum)[,-1])
        ## weights will change in oecosimu thus need to be recalculated
        if (weights == "prop")
            wt <- lapply(1:nl, function(i) apply(tmp[[i]], 1, function(z) sum(z) / sum(matr)))
            else wt <- lapply(1:nl, function(i) rep(1 / NROW(tmp[[i]]), NROW(tmp[[i]])))
        a <- sapply(1:nl, function(i) sum(divfun(tmp[[i]]) * wt[[i]]))
        names(a) <- c(paste("alpha", 1:(nl-1), sep="."), "gamma")
        b <- sapply(2:nl, function(i) a[i] - a[(i-1)])
        names(b) <- paste("beta", 1:(nl-1), sep=".")
        c(a, b)
    }
    sim <- oecosimu(matr, wdivfun, method = "permat", nsimul=nsimul,
        burnin=control$burnin, thin=control$thin, control=control)
    attr(sim, "index") <- index
    attr(sim, "weights") <- weights
    attr(sim, "n.levels") <- nl
    class(sim) <- c("adipart", "list")
    sim
}
