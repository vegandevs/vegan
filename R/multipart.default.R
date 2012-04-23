multipart.default <-
function(y, x, index=c("renyi", "tsallis"), scales = 1,
    global = FALSE, relative = FALSE, nsimul=99, ...)
{
    if (length(scales) > 1)
        stop("length of 'scales' must be 1")
    ## evaluate formula
    lhs <- as.matrix(y)
    if (missing(x))
        x <- cbind(level_1=seq_len(nrow(lhs)), 
            leve_2=rep(1, nrow(lhs)))
    rhs <- data.frame(x)
    rhs[] <- lapply(rhs, as.factor)
    rhs[] <- lapply(rhs, droplevels)
    nlevs <- ncol(rhs)
    if (nlevs < 2)
        stop("provide at least two level hierarchy")
    if (any(rowSums(lhs) == 0))
        stop("data matrix contains empty rows")
    if (any(lhs < 0))
        stop("data matrix contains negative entries")
    if (is.null(colnames(rhs)))
        colnames(rhs) <- paste("level", 1:nlevs, sep="_")
    tlab <- colnames(rhs)

     ## part check proper design of the model frame
    l1 <- sapply(rhs, function(z) length(unique(z)))
    if (!any(sapply(2:nlevs, function(z) l1[z] <= l1[z-1])))
        stop("number of levels are inapropriate, check sequence")
    rval <- list()
    rval[[1]] <- rhs[,nlevs]
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
#    if (!fullgamma && !global)
#        warning("gamma diversity value might be meaningless")
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
    index <- match.arg(index)
    divfun <- switch(index,
        "renyi" = function(x) renyi(x, scales=scales, hill = TRUE),
        "tsallis" = function(x) tsallis(x, scales=scales, hill = TRUE))

    ## cluster membership determination
    nrhs <- rhs
    nrhs <- sapply(nrhs, as.numeric)
    idcl <- function(i) {
        h <- nrhs[,i]
        l <- nrhs[,(i-1)]
        sapply(unique(l), function(i) h[l==i][1])
    }
    id <- lapply(2:nlevs, idcl)

    ## this is the function passed to oecosimu
    if (global) {
        wdivfun <- function(x) {
            if (fullgamma) {
                tmp <- lapply(1:(nlevs-1), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
                tmp[[nlevs]] <- matrix(colSums(x), nrow = 1, ncol = ncol(x))
            } else {
                tmp <- lapply(1:nlevs, function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
            }
            raw <- sapply(1:nlevs, function(i) divfun(tmp[[i]]))
            a <- sapply(raw, mean)
            G <- a[nlevs]
            b <- sapply(1:(nlevs-1), function(i) G / a[i])
            if (relative)
                b <- b / sapply(raw[1:(nlevs-1)], length)
            c(a, b)
        }
    } else {
        wdivfun <- function(x) {
            if (fullgamma) {
                tmp <- lapply(1:(nlevs-1), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
                tmp[[nlevs]] <- matrix(colSums(x), nrow = 1, ncol = ncol(x))
            } else {
                tmp <- lapply(1:nlevs, function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
            }
            a <- sapply(1:nlevs, function(i) divfun(tmp[[i]]))
            am <- lapply(1:(nlevs-1), function(i) {
                    sapply(1:length(unique(id[[i]])), function(ii) {
                        mean(a[[i]][id[[i]]==ii])
                    })
                })
            b <- lapply(1:(nlevs-1), function(i) a[[(i+1)]] / am[[i]])
            bmax <- lapply(id, function(i) table(i))
            if (relative)
                b <- lapply(1:(nlevs-1), function(i) b[[i]] / bmax[[i]])
            c(sapply(a, mean), sapply(b, mean))
        }
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
    attr(sim, "scales") <- scales
    attr(sim, "global") <- TRUE
    attr(sim, "n.levels") <- nlevs
    attr(sim, "terms") <- tlab
    attr(sim, "model") <- rhs
    class(sim) <- c("multipart", "list")
    sim
}
