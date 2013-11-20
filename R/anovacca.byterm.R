### Implementation of by-cases for vegan 2.2 versions of
### anova.cca. These are all internal functions that are not intended
### to be called by users in normal sessions, but they should be
### called from anova.cca (2.2). Therefore the user interface is rigid
### and input is not checked. The 'permutations' should be a
### permutation matrix.

### by = terms builds models as a sequence of adding terms and submits
### this to anova.ccalist

`anovacca.byterm` <-
    function(object, permutations, model, parallel)
{
    ## We need term labels but without Condition() terms
    trms <- terms(object)
    trmlab <- attr(trms, "term.labels")
    trmlab <- trmlab[trmlab %in% attr(terms(object$terminfo),
                                      "term.labels")]
    m0 <- update(object, . ~ 1)
    mods <- list(m0)
    for(i in seq_along(trmlab)) {
        fla <- paste(". ~ . + ", trmlab[i])
        mods[[i+1]] <- update(mods[[i]], fla)
    }
    ## The result. Should be reformatted
    anova.ccalist(mods, permutations = permutations, model = model, parallel = parallel)
}

## by = margin: this is not a anova.ccalist case, but we omit each
## term in turn and compare against the complete model.

`anovacca.bymargin` <-
    function(object, permutations, ...)
{
    nperm <- nrow(permutations)
    ## We need term labels but without Condition() terms
    trms <- terms(object)
    trmlab <- attr(trms, "term.labels")
    trmlab <- trmlab[trmlab %in% attr(terms(object$terminfo),
                                      "term.labels")]
    ## baseline: all terms
    big <- permutest(object, permutations, ...)
    dfbig <- big$df[2]
    chibig <- big$chi[2]
    scale <- big$den/dfbig
    ## Collect all marginal models
    mods <- lapply(trmlab, function(nm, ...)
           permutest(update(object, paste(".~.-", nm)),
                     permutations, ...))
    ## Chande in df
    Df <- sapply(mods, function(x) x$df[2]) - dfbig
    ## F of change
    Chisq <- sapply(mods, function(x) x$chi[2]) - chibig
    Fstat <- (Chisq/Df)/(chibig/dfbig)
    ## Simulated F-values
    Fval <- sapply(mods, function(x) x$den)
    Fval <- sweep(Fval, 1, big$den, "-")
    Fval <- sweep(Fval, 2, Df, "/")
    Fval <- sweep(Fval, 1, scale, "/")
    ## Simulated P-values
    Pval <- (colSums(sweep(Fval, 2, Fstat, ">=")) + 1)/(nperm + 1)
    ## Collect results to anova data.frame
    out <- data.frame(c(Df, dfbig), c(Chisq, chibig),
                      c(Fstat, NA), c(Pval, NA))
    colnames(out) <- c("Df", "Chisq", "F", "Pr(>F)")
    rownames(out) <- c(trmlab, "Residual")
    class(out) <- c("anova", "data.frame")
    out
}
