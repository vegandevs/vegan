### Implementation of by-cases for anova.cca. These are all internal
### functions that are not intended to be called by users in normal
### sessions, but they should be called from anova.cca. Therefore the
### user interface is rigid and input is not checked. The
### 'permutations' should be a permutation matrix.

### by = "terms" calls directly permutest.cca which decomposes the
### inertia between successive terms within compiled C code.

`anovaCCAbyterm` <-
    function(object, permutations, model, parallel)
{
    ## The result
    sol <- permutest(object, permutations = permutations,
                     model = model, by = "terms", parallel = parallel)
    ## Reformat
    EPS <- sqrt(.Machine$double.eps)
    Pval <- (colSums(sweep(sol$F.perm, 2, sol$F.0 - EPS, ">=")) + 1) /
        (sol$nperm + 1)
    out <- data.frame(sol$df, sol$chi, c(sol$F.0, NA), c(Pval, NA))

    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    dimnames(out) <- list(c(sol$termlabels, "Residual"),
                          c("Df", varname, "F", "Pr(>F)"))
    head <- paste0("Permutation test for ", object$method, " under ",
                   model, " model\n",
                   "Terms added sequentially (first to last)\n",
                   howHead(attr(permutations, "control")))
    mod <- paste("Model:", c(object$call))
    attr(out, "heading") <- c(head, mod)
    attr(out, "F.perm") <- sol$F.perm
    class(out) <- c("anova.cca", "anova","data.frame")
    out
}

## by = "margin": we omit each term in turn and compare against the
## complete model. This does not involve partial terms (Conditions) on
## other variables, but the permutations remain similar in "direct"
## and "reduced" (default) models (perhaps this model should not be
## used with "full" models?). This is basically similar decomposition
## as by="term", but compares models without each term in turn against
## the complete model in separate calls to permutest.cca. From vegan
## 2.5-0 this does not update model formula -- this avoids scoping
## issues and makes the function more robust when embedded in other
## functions. Instead, we call ordConstrained with method="pass" with
## modified constraint matrix.

`anovaCCAbymargin` <-
    function(object, permutations, scope, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    nperm <- nrow(permutations)
    ## We need term labels but without Condition() terms
    if (!is.null(scope) && is.character(scope))
        trms <- scope
    else
        trms <- drop.scope(object)
    ## Condition() not considered marginal
    alltrms <- intersect(attr(terms(object$terminfo), "term.labels"),
                         attr(terms(object), "term.labels"))
    trmlab <- intersect(alltrms, trms)
    if (length(trmlab) == 0)
        stop("the scope was empty: no available marginal terms")
    ## baseline: all terms
    big <- permutest(object, permutations, ...)
    dfbig <- big$df[2]
    chibig <- big$chi[2]
    scale <- big$den/dfbig
    ## Collect all marginal models. This differs from old version
    ## (vegan 2.0) where other but 'nm' were partialled out within
    ## Condition(). Now we only fit the model without 'nm' and compare
    ## the difference against the complete model.
    Y <- ordiYbar(object, "init")
    X <- model.matrix(object)
    ## we must have Constraints to get here, but we may also have
    ## Conditions
    if (!is.null(object$pCCA)) {
        Z <- X$Conditions
        X <- X$Constraints
    } else {
        Z <- NULL
    }
    ass <- object$terminfo$assign
    if (is.null(ass))
        stop("old style result object: update() your model")
    ## analyse only terms of 'ass' thar are in scope
    scopeterms <- which(alltrms %in% trmlab)
    mods <- lapply(scopeterms, function(i, ...)
           permutest(ordConstrained(Y, X[, ass != i, drop=FALSE], Z, "pass"),
                     permutations, ...), ...)
    ## Chande in df
    Df <- sapply(mods, function(x) x$df[2]) - dfbig
    ## F of change
    Chisq <- sapply(mods, function(x) x$chi[2]) - chibig
    Fstat <- (Chisq/Df)/(chibig/dfbig)
    ## Simulated F-values
    Fval <- sapply(mods, function(x) x$num)
    ## Had we an empty model we need to clone the denominator
    if (length(Fval) == 1)
        Fval <- matrix(Fval, nrow = nperm)
    Fval <- sweep(-Fval, 1, big$num, "+")
    Fval <- sweep(Fval, 2, Df, "/")
    Fval <- sweep(Fval, 1, scale, "/")
    ## Simulated P-values
    Pval <- (colSums(sweep(Fval, 2, Fstat - EPS, ">=")) + 1)/(nperm + 1)
    ## Collect results to anova data.frame
    out <- data.frame(c(Df, dfbig), c(Chisq, chibig),
                      c(Fstat, NA), c(Pval, NA))
    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    dimnames(out) <- list(c(trmlab, "Residual"),
                          c("Df", varname, "F", "Pr(>F)"))
    head <- paste0("Permutation test for ", object$method, " under ",
                   big$model, " model\n",
                   "Marginal effects of terms\n",
                   howHead(attr(permutations, "control")))
    mod <- paste("Model:", c(object$call))
    attr(out, "heading") <- c(head, mod)
    attr(out, "F.perm") <- Fval
    class(out) <- c("anova.cca", "anova", "data.frame")
    out
}

### by = "axis" uses partial model: we use the original constraints,
### but add previous axes 1..(k-1) to Conditions when evaluating the
### significance of axis k which is compared against the first
### eigenvalue of the permutations. To avoid scoping issues, this
### calls directly ordConstrained() with modified Conditions (Z) and
### original Constraints (X) instead of updating formula. This
### corresponds to "forward" model in Legendre, Oksanen, ter Braak
### (2011).

### In 2.2-x to 2.4-3 we used "marginal model" where original
### Constraints were replaced with LC scores axes (object$CCA$u), and
### all but axis k were used as Conditions when evaluating the
### significance of axis k. My (J.Oksanen) simulations showed that
### this gave somewhat biased results.

`anovaCCAbyaxis` <-
    function(object, permutations, model, parallel, cutoff = 1)
{
    EPS <- sqrt(.Machine$double.eps)
    ## On 29/10/15 (983ba7726) we assumed that dbrda(d ~ dbrda(d ~
    ## x)$CCA$u) is not equal to dbrda(d ~ x) when there are negative
    ## eigenvalues, but it seems that it is OK if constrained
    ## eigenvalues are non-negative
    if (inherits(object, "dbrda") && any(object$CCA$eig < 0))
        stop("by = 'axis' cannot be used when constraints have negative eigenvalues")
    nperm <- nrow(permutations)
    ## Observed F-values and Df
    eig <- object$CCA$eig
    resdf <- nobs(object) - length(eig) - max(object$pCCA$QR$rank, 0) - 1
    Fstat <- eig/object$CA$tot.chi*resdf
    Df <- rep(1, length(eig))

    ## collect header and varname here: 'object' is modified later
    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"

    head <- paste0("Permutation test for ", object$method, " under ",
                   model, " model\n",
                   "Forward tests for axes\n",
                   howHead(attr(permutations, "control")))
    head <- c(head, paste("Model:", c(object$call)))

    ## constraints and model matrices
    Y <- object$Ybar
    if (is.null(Y))
        stop("old style result object: update() your model")
    if (!is.null(object$pCCA))
        Z <- qr.X(object$pCCA$QR)
    else
        Z <- NULL
    X <- model.matrix(object)
    if (!is.null(object$pCCA)) {
        Z <- X$Conditions
        X <- X$Constraints
    } else {
        Z <- NULL
    }
    LC <- object$CCA$u

    Pvals <- rep(NA, ncol(LC))
    F.perm <- matrix(ncol = ncol(LC), nrow = nperm)
    for (i in seq_along(eig)) {
        if (i > 1) {
            object <- ordConstrained(Y, X, cbind(Z, LC[, seq_len(i-1)]), "pass")
        }
        if (length(eig) == i) {
            mod <- permutest(object, permutations, model = model,
                             parallel = parallel)
        } else {
            mod <- permutest(object, permutations, model = model,
                             parallel = parallel, first = TRUE)
        }
        Pvals[i] <- (sum(mod$F.perm >= mod$F.0 - EPS) + 1) / (nperm + 1)
        ## follow Canoco: P-values of later axes cannot be lower than
        ## previous axes (usually no effect as P-values are increasing).
        if (i > 1 && Pvals[i] < Pvals[i-1])
            Pvals[i] <- Pvals[i-1]
        F.perm[ , i] <- mod$F.perm
        if (Pvals[i] >= cutoff)
            break
    }
    out <- data.frame(c(Df, resdf), c(eig, object$CA$tot.chi),
                      c(Fstat, NA), c(Pvals,NA))
    rownames(out) <- c(names(eig), "Residual")
    colnames(out) <- c("Df", varname, "F", "Pr(>F)")
    attr(out, "heading") <- head
    attr(out, "F.perm") <- F.perm
    class(out) <- c("anova.cca", "anova", "data.frame")
    out
}

### Wrap permutest.cca(..., by="onedf") in a anova.cca form

`anovaCCAby1df` <-
    function(object, permutations, model, parallel)
{
    ## Compute
    sol <- permutest(object, permutations = permutations,
                     model = model, by = "onedf", parallel = parallel)
    ## Reformat
    EPS <- sqrt(.Machine$double.eps)
    Pval <- (colSums(sweep(sol$F.perm, 2, sol$F.0 - EPS, ">=")) + 1) /
        (sol$nperm + 1)
    out <- data.frame(sol$df, sol$chi, c(sol$F.0, NA), c(Pval, NA))

    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    dimnames(out) <- list(c(sol$termlabels, "Residual"),
                          c("Df", varname, "F", "Pr(>F)"))
    head <- paste0("Permutation test for ", object$method, " under ",
                   model, " model\n",
                   "Sequential test for contrasts\n",
                   howHead(attr(permutations, "control")))
    mod <- paste("Model:", c(object$call))
    attr(out, "heading") <- c(head, mod)
    attr(out, "F.perm") <- sol$F.perm
    class(out) <- c("anova.cca", "anova","data.frame")
    out
}
