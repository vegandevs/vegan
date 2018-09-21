### The constrained ordination methods cca, rda, capscale & dbrda are
### very similar. Their main difference is the init of dependent
### matrix and after that the partializing, constraining and residual
### steps are very similar to each other. This file provides functions
### that can analyse any classic method, the only difference being the
### attributes of the dependet variable.

### In this file we use convention modelfunction(Y, X, Z) where Y is
### the dependent data set (community), X is the model matrix of
### constraints, and Z is the model matrix of conditions to be
### partialled out. The common outline of the function is:
###
### initX(Y)
### pCCA <- ordPartial(Y, Z)
### CCA <- ordConstrain(Y, X, Z)
### CA <- ordResid(Y)
###
### The init step sets up the dependent data Y and sets up its
### attributes that control how it will be handled at later
### stage. Each of the later stages modifies Y and returns results of
### its analysis. The handling of the function is mainly similar, but
### there is some variation depending on the attributes of Y.

### THE USAGE

### Function prototype is ordConstrained(Y, X=NULL, Z=NULL, method),
### where Y is the dependent community data set, X is the model matrix
### of constraints, Z the model matrix of conditions, method is "cca",
### "rda", "capscale" or "dbrda" (the last two may not work, and
### "capscale" may not work in the way you assume). The function
### returns a large subset of correspoding constrained ordination
### method. For instance, with method = "dbrda", the result is mainly
### correct, but it differs so much from the current dbrda that it
### cannot be printed cleanly.

### THE INIT METHODS

### The init methods transform the dependent data specifically to the
### particular method and set up attributes that will control further
### processing of Y. The process critical attributes are set up in
### UPPER CASE to make the if-statements stand out in later analysis.

`initPCA` <-
    function(Y, scale = FALSE)
{
    Y <- as.matrix(Y)
    Y <- scale(Y, scale = scale)
    if (scale && any(is.nan(Y)))
        Y[is.nan(Y)] <- 0
    ## we want models based on variance or correlations -- this will
    ## break methods depending on unscaled Xbar (and which do this
    ## very same scaling internally) with scale = FALSE.
    Y <- Y / sqrt(nrow(Y) - 1)
    attr(Y, "METHOD") <- "PCA"
    Y
}

`initCAP` <-
    function(Y)
{
    Y <- as.matrix(Y)
    Y <- scale(Y, scale = FALSE)
    attr(Y, "METHOD") <- "CAPSCALE"
    Y
}

`initCA` <-
    function(Y)
{
    Y <- as.matrix(Y)
    tot <- sum(Y)
    Y <- Y/tot
    rw <- rowSums(Y)
    cw <- colSums(Y)
    rc <- outer(rw, cw)
    Y <- (Y - rc)/sqrt(rc)
    attr(Y, "tot") <- tot
    attr(Y, "RW") <- rw
    attr(Y, "CW") <- cw
    attr(Y, "METHOD") <- "CA"
    Y
}

`initDBRDA` <-
    function(Y)
{
    ## check
    Y <- as.matrix(Y)
    dims <- dim(Y)
    if (dims[1] != dims[2] || !isSymmetric(unname(Y)))
        stop("input Y must be distances or a symmetric square matrix")
    ## transform
    Y <- -0.5 * GowerDblcen(Y^2)
    attr(Y, "METHOD") <- "DISTBASED"
    Y
}

### COMMON HEADER INFORMATION FOR ORDINATION MODELS

`ordHead`<- function(Y)
{
    method <- attr(Y, "METHOD")
    headmethod <- switch(method,
                         "CA" = "cca",
                         "PCA" = "rda",
                         "CAPSCALE" = "capscale",
                         "DISTBASED" = "dbrda")
    if (method == "DISTBASED")
        totvar <- sum(diag(Y))
    else
        totvar <- sum(Y^2)
    head <- list("tot.chi" = totvar, "Ybar" = Y, "method" = headmethod)
    if (method == "CA")
        head <- c(list("grand.total" = attr(Y, "tot"),
                       "rowsum" = attr(Y, "RW"),
                       "colsum" = attr(Y, "CW")),
                  head)
    else if (method == "PCA")
        head <- c(list("colsum" = sqrt(colSums(Y^2))),
                  head)
    head
}

### THE PARTIAL MODEL

`ordPartial` <-
    function(Y, Z)
{
    ZERO <- sqrt(.Machine$double.eps)
    ## attributes
    DISTBASED <- attr(Y, "METHOD") == "DISTBASED"
    RW <- attr(Y, "RW")
    ## centre Z
    if (!is.null(RW)) {
        envcentre <- apply(Z, 2, weighted.mean, w = RW)
        Z <- scale(Z, center = envcentre, scale = FALSE)
        Z <- sweep(Z, 1, sqrt(RW), "*")
    } else {
        envcentre <- colMeans(Z)
        Z <- scale(Z, center = envcentre, scale = FALSE)
    }
    ## QR decomposition
    Q <- qr(Z)
    ## partialled out variation as a trace of Yfit
    Yfit <- qr.fitted(Q, Y)
    if (DISTBASED) {
        Yfit <- qr.fitted(Q, t(Yfit))
        totvar <- sum(diag(Yfit))
    } else {
        totvar <- sum(Yfit^2)
    }
    if (totvar < ZERO)
        totvar <- 0
    ## residuals of Y
    Y <- qr.resid(Q, Y)
    if (DISTBASED)
        Y <- qr.resid(Q, t(Y))
    ## result object like in current cca, rda
    result <- list(
        rank = if (totvar > 0) Q$rank else 0,
        tot.chi = totvar,
        QR = Q,
        Fit = Yfit,
        envcentre = envcentre)
    list(Y = Y, result = result)
}

### THE CONSTRAINTS

`ordConstrain` <- function(Y, X, Z)
{
    ## attributes & constants
    DISTBASED <- attr(Y, "METHOD") == "DISTBASED"
    RW <- attr(Y, "RW")
    CW <- attr(Y, "CW")
    ZERO <- sqrt(.Machine$double.eps)
    ## combine conditions and constraints if necessary
    if (!is.null(Z)) {
        X <- cbind(Z, X)
        zcol <- ncol(Z)
    } else {
        zcol <- 0
    }
    ## centre
    if (!is.null(RW)) {
        envcentre <- apply(X, 2, weighted.mean, w = RW)
        X <- scale(X, center = envcentre, scale = FALSE)
        X <- sweep(X, 1, sqrt(RW), "*")
    } else {
        envcentre <- colMeans(X)
        X <- scale(X, center = envcentre, scale = FALSE)
    }
    ## QR
    Q <- qr(X)
    ## we need to see how much rank grows over rank of conditions
    rank <- sum(Q$pivot[seq_len(Q$rank)] > zcol)
    ## check for aliased terms
    if (length(Q$pivot) > Q$rank)
        alias <- colnames(Q$qr)[-seq_len(Q$rank)]
    else
        alias <- NULL
    ## kept constraints and their means
    kept <- seq_along(Q$pivot) <= Q$rank & Q$pivot > zcol
    if (zcol > 0)
        envcentre <- envcentre[-seq_len(zcol)]
    ## eigen solution
    Yfit <- qr.fitted(Q, Y)
    if (DISTBASED) {
        Yfit <- qr.fitted(Q, t(Yfit))
        sol <- eigen(Yfit, symmetric = TRUE)
        lambda <- sol$values
        u <- sol$vectors
    } else {
        sol <- svd(Yfit)
        lambda <- sol$d^2
        u <- sol$u
        v <- sol$v
    }
    ## handle zero  eigenvalues and negative eigenvalues
    zeroev <- abs(lambda) < max(ZERO, ZERO * lambda[1L])
    if (any(zeroev)) {
        lambda <- lambda[!zeroev]
        u <- u[, !zeroev, drop = FALSE]
        if (!DISTBASED)
            v <- v[, !zeroev, drop = FALSE]
    }
    posev <- lambda > 0
    ## wa scores
    if (DISTBASED) {
        wa <- Y %*% u[, posev, drop = FALSE] %*%
            diag(1/lambda[posev], sum(posev))
        v <- matrix(NA, 0, sum(posev))
    } else {
        wa <- Y %*% v %*% diag(1/sqrt(lambda), sum(posev))
    }
    ## biplot scores: basically these are cor(X, u), but cor() would
    ## re-centre X and u in CCA, and therefore we need the following
    ## (which also is faster)
    xx <- X[, Q$pivot[kept], drop = FALSE]
    bp <- (1/sqrt(colSums(xx^2))) * crossprod(xx, u[, posev, drop=FALSE])
    ## de-weight
    if (!is.null(RW)) {
        u <- sweep(u, 1, sqrt(RW), "/")
        if (all(!is.na(wa)))
            wa <- sweep(wa, 1, sqrt(RW), "/")
    }
    if (!is.null(CW) && nrow(v)) {
        v <- sweep(v, 1, sqrt(CW), "/")
    }
    ## set names
    axnam <- paste0(switch(attr(Y, "METHOD"),
                           "PCA" = "RDA",
                           "CA" = "CCA",
                           "CAPSCALE" = "CAP",
                           "DISTBASED" = "dbRDA"),
                    seq_len(sum(posev)))
    if (DISTBASED && any(!posev))
        negnam <- paste0("idbRDA", seq_len(sum(!posev)))
    else
        negnam <- NULL
    dnam <- dimnames(Y)
    if (any(posev))
        names(lambda) <- c(axnam, negnam)
    if (ncol(u))
        dimnames(u) <- list(dnam[[1]], c(axnam, negnam))
    if (nrow(v) && ncol(v))     # no rows in DISTBASED
        dimnames(v) <- list(dnam[[2]], axnam)
    if (ncol(wa)) # only for posev
        colnames(wa) <- axnam
    if (ncol(bp))    # only for posev
        colnames(bp) <- axnam
    ## out
    result <- list(
        eig = lambda,
        poseig = if (DISTBASED) sum(posev) else NULL,
        u = u,
        v = v,
        wa = wa,
        alias = alias,
        biplot = bp,
        rank = length(lambda),
        qrank = rank,
        tot.chi = sum(lambda),
        QR = Q,
        envcentre = envcentre,
        Xbar = Y)
    ## residual of Y
    Y <- qr.resid(Q, Y)
    if (DISTBASED)
        Y <- qr.resid(Q, t(Y))
    ## out
    list(Y = Y, result = result)
}

### THE RESIDUAL METHOD

### Finds the unconstrained ordination after (optionally) removing the
### variation that could be explained by partial and constrained
### models.

`ordResid` <-
    function(Y)
{
    ## get attributes
    DISTBASED <- attr(Y, "METHOD") == "DISTBASED"
    RW <- attr(Y, "RW")
    CW <- attr(Y, "CW")
    ## Ordination
    ZERO <- sqrt(.Machine$double.eps)
    if (DISTBASED) {
        sol <- eigen(Y, symmetric = TRUE)
        lambda <- sol$values
        u <- sol$vectors
    } else {
        sol <- svd(Y)
        lambda <- sol$d^2
        u <- sol$u
        v <- sol$v
    }
    ## handle zero and negative eigenvalues
    zeroev <- abs(lambda) < max(ZERO, ZERO * lambda[1L])
    if (any(zeroev)) {
        lambda <- lambda[!zeroev]
        u <- u[, !zeroev, drop = FALSE]
        if (!DISTBASED) # no v in DISTBASED
            v <- v[, !zeroev, drop = FALSE]
    }
    posev <- lambda > 0
    if (DISTBASED) # no species scores in DISTBASED
        v <- matrix(NA, 0, sum(posev))
    ## de-weight
    if (!is.null(RW)) {
        u <- sweep(u, 1, sqrt(RW), "/")
    }
    if (!is.null(CW) && nrow(v)) {
        v <- sweep(v, 1, sqrt(CW), "/")
    }
    ## set names
    axnam <- paste0(switch(attr(Y, "METHOD"),
                           "PCA" = "PC",
                           "CA" = "CA",
                           "CAPSCALE" = "MDS",
                           "DISTBASED" = "MDS"),
                    seq_len(sum(posev)))
    if (DISTBASED && any(!posev))
        negnam <- paste0("iMDS", seq_len(sum(!posev)))
    else
        negnam <- NULL
    dnam <- dimnames(Y)
    if (any(posev))
        names(lambda) <- c(axnam, negnam)
    if (ncol(u))
        dimnames(u) <- list(dnam[[1]], c(axnam, negnam))
    if (nrow(v) && ncol(v)) # no rows in DISTBASED
        dimnames(v) <- list(dnam[[2]], axnam)
    ## out
    out <- list(
        "eig" = lambda,
        "poseig" = if (DISTBASED) sum(posev) else NULL,
        "u" = u,
        "v" = v,
        "rank" = length(lambda),
        "tot.chi" = sum(lambda),
        "Xbar" = Y)
    out
}

## The actual function that calls all previous and returns the fitted
## ordination model

`ordConstrained` <-
    function(Y, X = NULL, Z = NULL,
             method = c("cca", "rda", "capscale", "dbrda", "pass"),
             arg = FALSE)
{
    method = match.arg(method)
    partial <- constraint <- resid <- NULL
    ## init; "pass" returns unchanged Y, presumably from previous init
    Y <- switch(method,
                "cca" = initCA(Y),
                "rda" = initPCA(Y, scale = arg),
                "capscale" = initCAP(Y),
                "dbrda" = initDBRDA(Y),
                "pass" = Y)
    ## header info for the model
    head <- ordHead(Y)
    ## Partial
    if (!is.null(Z) && ncol(Z)) {
        out <- ordPartial(Y, Z)
        Y <- out$Y
        partial <- out$result
    }
    ## Constraints
    if (!is.null(X) && ncol(X)) {
        out <- ordConstrain(Y, X, Z)
        Y <- out$Y
        constraint <- out$result
    }
    ## Residuals
    resid <- ordResid(Y)
    ## return a CCA object
    out <- c(head,
             call = match.call(),
             list("pCCA" = partial, "CCA" = constraint, "CA" = resid))
    class(out) <- switch(attr(Y, "METHOD"),
                         "CA" = "cca",
                         "PCA" = c("rda", "cca"),
                         "CAPSCALE" = c("capscale", "rda", "cca"),
                         "DISTBASED" = c("dbrda", "rda", "cca"))
    out
}
