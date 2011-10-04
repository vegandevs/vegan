permutest <- function(x, ...)
    UseMethod("permutest")

permutest.default <- function(x, ...)
    stop("No default permutation test defined")

`permutest.cca` <-
    function (x, permutations = 99,
              model = c("reduced", "direct", "full"), first = FALSE,
              strata, ...) 
{
    model <- match.arg(model)
    isCCA <- !inherits(x, "rda")
    isPartial <- !is.null(x$pCCA)
    if (first) {
        Chi.z <- x$CCA$eig[1]
        q <- 1
    }
    else {
        Chi.z <- x$CCA$tot.chi
        names(Chi.z) <- "Model"
        q <- x$CCA$qrank
    }
    Chi.xz <- x$CA$tot.chi
    names(Chi.xz) <- "Residual"
    r <- nrow(x$CA$Xbar) - x$CCA$QR$rank - 1
    if (model == "full") 
        Chi.tot <- Chi.xz
    else Chi.tot <- Chi.z + Chi.xz
    if (!isCCA) 
        Chi.tot <- Chi.tot * (nrow(x$CCA$Xbar) - 1)
    F.0 <- (Chi.z/q)/(Chi.xz/r)
    F.perm <- numeric(permutations)
    num <- numeric(permutations)
    den <- numeric(permutations)
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
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    for (i in 1:permutations) {
        take <- permuted.index(N, strata)
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
        num[i] <- cca.ev
        den[i] <- ca.ev
        F.perm[i] <- (cca.ev/q)/(ca.ev/r)
    }
    ## Round to avoid arbitrary ordering of statistics due to
    ## numerical inaccuracy
    F.0 <- round(F.0, 12)
    F.perm <- round(F.perm, 12)
    Call <- match.call()
    Call[[1]] <- as.name("permutest")
    sol <- list(call = Call, testcall = x$call, model = model,
                F.0 = F.0, F.perm = F.perm,  chi = c(Chi.z, Chi.xz),
                num = num, den = den, df = c(q, r), nperm = permutations,
                method = x$method, first = first,  Random.seed = seed)
    if (!missing(strata)) {
        sol$strata <- deparse(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "permutest.cca"
    sol
}
