### Implementation of by-cases for vegan 2.2 versions of
### anova.cca. These are all internal functions that are not intended
### to be called by users in normal sessions, but they should be
### called from anova.cca. Therefore the user interface is rigid and
### input is not checked. The 'permutations' should be a permutation
### matrix.

### by = terms calls directly permutest.cca

`anova.ccabyterm` <-
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

## by = margin: this is not a anova.ccalist case, but we omit each
## term in turn and compare against the complete model.

`anova.ccabymargin` <-
    function(object, permutations, scope, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    nperm <- nrow(permutations)
    ## Refuse to handle models with missing data
    if (!is.null(object$na.action))
        stop("by = 'margin' models cannot handle missing data")
    ## We need term labels but without Condition() terms
    if (!is.null(scope) && is.character(scope))
        trms <- scope
    else
        trms <- drop.scope(object)
    trmlab <- trms[trms %in% attr(terms(object$terminfo),
                                      "term.labels")]
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
    mods <- lapply(trmlab, function(nm, ...)
           permutest(update(object, paste(".~.-", nm)),
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
                   mods[[1]]$model, " model\n",
                   "Marginal effects of terms\n",
                   howHead(attr(permutations, "control")))
    mod <- paste("Model:", c(object$call))
    attr(out, "heading") <- c(head, mod)
    attr(out, "F.perm") <- Fval
    class(out) <- c("anova.cca", "anova", "data.frame")
    out
}

### Marginal test for axes

`anova.ccabyaxis` <-
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
    resdf <- nobs(object) - length(eig) - max(object$pCCA$rank, 0) - 1
    Fstat <- eig/object$CA$tot.chi*resdf
    Df <- rep(1, length(eig))
    ## save object: will be modified later
    origobj <- object
    ## constraints and model matrices
    Y <- object$Ybar
    if (is.null(Y))
        stop("old style result object does not work: update() your model")
    if (!is.null(object$pCCA))
        Z <- qr.X(object$pCCA$QR)
    else
        Z <- NULL
    X <- qr.X(object$CCA$QR)
    LC <- object$CCA$u
    ## In CA we need to de-weight X and Z
    if (attr(Y, "METHOD") == "CA") {
        invw <- 1/sqrt(attr(Y, "RW"))
        if (!is.null(Z))
            Z <- invw * Z
        X <- invw * X
    }
    Pvals <- rep(NA, ncol(LC))
    F.perm <- matrix(ncol = ncol(LC), nrow = nperm)
    axnams <- colnames(LC)
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
        F.perm[ , i] <- mod$F.perm
        if (Pvals[i] > cutoff)
            break
    }
    out <- data.frame(c(Df, resdf), c(eig, origobj$CA$tot.chi),
                      c(Fstat, NA), c(Pvals,NA))
    rownames(out) <- c(names(eig), "Residual")
    if (inherits(origobj, c("capscale", "dbrda")) && origobj$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(origobj, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    colnames(out) <- c("Df", varname, "F", "Pr(>F)")
    head <- paste0("Permutation test for ", origobj$method, " under ",
                   model, " model\n",
                   "Marginal tests for axes\n",
                   howHead(attr(permutations, "control")))
    mod <- paste("Model:", c(origobj$call))
    attr(out, "heading") <- c(head, mod)
    attr(out, "F.perm") <- F.perm
    class(out) <- c("anova.cca", "anova", "data.frame")
    out
}
