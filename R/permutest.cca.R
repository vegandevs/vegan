permutest <- function(x, ...)
    UseMethod("permutest")

permutest.default <- function(x, ...)
    stop("No default permutation test defined")

`permutest.cca` <-
    function (x, permutations = how(nperm=99),
              model = c("reduced", "direct"), first = FALSE,
              strata = NULL, parallel = getOption("mc.cores") , ...)
{
    ## do something sensible with insensible input (no constraints)
    if (is.null(x$CCA)) {
        sol <- list(call = match.call(), testcall = x$call, model = NA,
                    F.0 = NA, F.perm = NA, chi = c(0, x$CA$tot.chi),
                    num = 0, den = x$CA$tot.chi,
                    df = c(0, nrow(x$CA$u) - max(x$pCCA$rank,0) - 1),
                    nperm = 0, method = x$method, first = FALSE,
                    Random.seed = NA)
        class(sol) <- "permutest.cca"
        return(sol)
    }
    model <- match.arg(model)
    isCCA <- !inherits(x, "rda")
    isPartial <- !is.null(x$pCCA)
    ## Function to get the F statistics in one loop
    getF <- function (indx, ...)
    {
        if (!is.matrix(indx))
            dim(indx) <- c(1, length(indx))
        R <- nrow(indx)
        mat <- matrix(0, nrow = R, ncol = 3)
        for (i in seq_len(R)) {
            take <- indx[i,]
            Y <- E[take, ]
            if (isCCA)
                wtake <- w[take]
            if (isPartial) {
                if (isCCA) {
                    XZ <- .C("wcentre", x = as.double(Z), as.double(wtake),
                             as.integer(N), as.integer(Zcol),
                             PACKAGE = "vegan")$x
                    dim(XZ) <- c(N, Zcol)
                    QZ <- qr(XZ)
                }
                Y <- qr.resid(QZ, Y)
            }
            if (isCCA) {
                XY <- .C("wcentre", x = as.double(X), as.double(wtake),
                         as.integer(N), as.integer(Xcol),
                         PACKAGE = "vegan")$x
                dim(XY) <- c(N, Xcol)
                Q <- qr(XY)
            }
            tmp <- qr.fitted(Q, Y)
            if (first)
                cca.ev <- La.svd(tmp, nv = 0, nu = 0)$d[1]^2
            else cca.ev <- sum(tmp * tmp)
            if (isPartial || first) {
                tmp <- qr.resid(Q, Y)
                ca.ev <- sum(tmp * tmp)
            }
            else ca.ev <- Chi.tot - cca.ev
            mat[i,] <- cbind(cca.ev, ca.ev, (cca.ev/q)/(ca.ev/r))
        }
        mat
    }
    ## end getF()
    if (first) {
        Chi.z <- x$CCA$eig[1]
        q <- 1
    }
    else {
        Chi.z <- x$CCA$tot.chi
        names(Chi.z) <- "Model"
        q <- x$CCA$qrank
    }
    ## Set up
    Chi.xz <- x$CA$tot.chi
    names(Chi.xz) <- "Residual"
    r <- nrow(x$CA$Xbar) - x$CCA$QR$rank - 1
    if (model == "full")
        Chi.tot <- Chi.xz
    else Chi.tot <- Chi.z + Chi.xz
    if (!isCCA)
        Chi.tot <- Chi.tot * (nrow(x$CCA$Xbar) - 1)
    F.0 <- (Chi.z/q)/(Chi.xz/r)
    Q <- x$CCA$QR
    if (isCCA) {
        w <- x$rowsum # works with any na.action, weights(x) won't
        X <- qr.X(Q, ncol=length(Q$pivot))
        X <- sweep(X, 1, sqrt(w), "/")
    }
    if (isPartial) {
        Y.Z <- x$pCCA$Fit
        QZ <- x$pCCA$QR
        if (isCCA) {
            Z <- qr.X(QZ)
            Z <- sweep(Z, 1, sqrt(w), "/")
        }
    }
    if (model == "reduced" || model == "direct")
        E <- x$CCA$Xbar
    else E <- x$CA$Xbar
    if (isPartial && model == "direct")
        E <- E + Y.Z
    ## Save dimensions
    N <- nrow(E)
    if (isCCA) {
        Xcol <- ncol(X)
        if (isPartial)
            Zcol <- ncol(Z)
    }
    permutations <- getPermuteMatrix(permutations, N, strata = strata)
    nperm <- nrow(permutations)
    ## Parallel processing (similar as in oecosimu)
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if ((hasClus || parallel > 1)  && require(parallel)) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            tmp <- do.call(rbind,
                           mclapply(1:nperm,
                                    function(i) getF(permutations[i,]),
                                    mc.cores = parallel))
        } else {
            ## if hasClus, do not set up and stop a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
            }
            tmp <- parRapply(parallel, permutations, function(i) getF(i))
            tmp <- matrix(tmp, ncol=3, byrow=TRUE)
            if (!hasClus)
                stopCluster(parallel)
        }
    } else {
        tmp <- getF(permutations)
    }
    num <- tmp[,1]
    den <- tmp[,2]
    F.perm <- tmp[,3]
    ## Round to avoid arbitrary ordering of statistics due to
    ## numerical inaccuracy
    F.0 <- round(F.0, 12)
    F.perm <- round(F.perm, 12)
    Call <- match.call()
    Call[[1]] <- as.name("permutest")
    sol <- list(call = Call, testcall = x$call, model = model,
                F.0 = F.0, F.perm = F.perm,  chi = c(Chi.z, Chi.xz),
                num = num, den = den, df = c(q, r), nperm = nperm,
                method = x$method, first = first)
    sol$Random.seed <- attr(permutations, "seed")
    sol$control <- attr(permutations, "control")
    if (!missing(strata)) {
        sol$strata <- deparse(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "permutest.cca"
    sol
}
