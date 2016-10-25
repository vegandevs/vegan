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
    ## we want variance based model when scale = FALSE -- this will
    ## break Xbar where we want to have back the original scaling
    if (!scale)
        Y <- Y / sqrt(nrow(Y) - 1)
    attr(Y, "METHOD") <- "PCA"
    Y
}

`initCA` <-
    function(Y)
{
    Y <- as.matrix(Y)
    Y <- Y/sum(Y)
    rw <- rowSums(Y)
    cw <- colSums(Y)
    rc <- outer(rw, cw)
    Y <- (Y - rc)/sqrt(rc)
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
    warning("ordination header not yet implemented: print etc will fail")
}

### THE PARTIAL MODEL

`ordPartial` <-
    function(Y, Z)
{
    ## attributes
    DISTBASED <- attr(Y, "METHOD") == "DISTBASED"
    RW <- attr(Y, "RW")
    CW <- attr(Y, "CW")
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
    if (DISTBASED)
        totvar <- sum(diag(Yfit))
    else
        totvar <- sum(Yfit^2)
    ## residuals of Y
    Y <- qr.resid(Q, Y)
    if (DISTBASED)
        Y <- qr.resid(Q, t(Y))
    ## result object like in current cca, rda
    result <- list(
        rank = Q$rank,
        tot.chi = totvar,
        QR = Q,
        Fit = Yfit,
        envcentre = envcentre)
    list(Y = Y, result = result)
}

### THE CONSTRAINTS

`ordConstraints` <- function(Y, X, Z)
{
    if (is.null(X))
        stop("Constrained models are not yet implemented")
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
    rw <- attr(Y, "RW")
    cw <- attr(Y, "CW")
    ## Ordination
    if (DISTBASED) {
        sol <- eigen(Y)
        lambda <- sol$values
        u <- sol$vectors
        v <- NULL
    } else {
        sol <- svd(Y)
        lambda <- sol$d^2
        u <- sol$u
        v <- sol$v
    }
    ## handle zero and negative eigenvalues... not yet implemented

    ## de-weight
    if (!is.null(rw)) {
        u <- sweep(u, 1, sqrt(rw), "/")
    }
    if (!is.null(cw)) {
        v <- sweep(v, 1, sqrt(cw), "/")
    }
    ## out
    out <- list(
        "eig" = lambda,
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
             method = c("cca", "rda", "capscale", "dbrda"),
             scale = FALSE)
{
    method = match.arg(method)
    partial <- constraint <- resid <- NULL
    ## init
    Y <- switch(method,
                "cca" = initCA(Y),
                "rda" = initPCA(Y, scale = scale),
                "capscale" = initPCA(Y, scale = FALSE),
                "dbrda" = initDBRDA(Y))
    ## header info for the model
    head <- ordHead(Y)
    ## Partial
    if (!is.null(Z)) {
        out <- ordPartial(Y, Z)
        Y <- out$Y
        partial <- out$result
    }
    ## Constraints
    if (!is.null(X)) {
        out <- ordConstrained(Y, X, Z)
        Y <- out$Y
        constraint <- out$result
    }
    ## Residuals
    resid <- ordResid(Y)
    ## return a CCA object
    out <- list("pCCA" = partial, "CCA" = constraint, "CA" = resid)
    class(out) <- switch(method,
                         "cca" = "cca",
                         "rda" = c("rda", "cca"),
                         "capscale" = c("capscale", "rda", "cca"),
                         "dbrda" = c("dbrda", "rda", "cca"))
    out
}
