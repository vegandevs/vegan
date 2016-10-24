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
    ## we want variance based model when scale = FALSE
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

`ordConstrained` <-
    function(Y, X = NULL, Z = NULL,
             method = c("cca", "rda", "capscale", "dbrda"),
             scale = FALSE)
{
    method = match.arg(method)
    ## init
    Y <- switch(method,
                "cca" = initCA(Y),
                "rda" = initPCA(Y, scale = scale),
                "capscale" = initPCA(Y, scale = FALSE),
                "dbrda" = initDBRDA(Y))
    ## Partial
    if (!is.null(Z))
        stop("Partial models are not yet implemented")
    ## Constraints
    if (!is.null(X))
        stop("Constrained models are not yet implemented")
    ## Residuals
    ## not yet implemented: so return only Y
    Y
}
