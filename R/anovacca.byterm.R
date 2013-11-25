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
    ntrm <- length(trmlab)
    m0 <- update(object, paste(".~.-", paste(trmlab, collapse="-")))
    mods <- list(m0)
    for(i in seq_along(trmlab)) {
        fla <- paste(". ~ . + ", trmlab[i])
        mods[[i+1]] <- update(mods[[i]], fla)
    }
    ## The result
    sol <- anova.ccalist(mods, permutations = permutations,
                         model = model, parallel = parallel)
    ## Reformat
    out <- data.frame(c(sol[-1,3], sol[ntrm+1,1]),
                      c(sol[-1,4], sol[ntrm+1,2]),
                      c(sol[-1,5], NA),
                      c(sol[-1,6], NA))
    colnames(out) <- c("Df", "Chisq", "F", "Pr(>F)")
    rownames(out) <- c(trmlab, "Residual")
    class(out) <- c("anova","data.frame")
    out
}

## by = margin: this is not a anova.ccalist case, but we omit each
## term in turn and compare against the complete model.

`anovacca.bymargin` <-
    function(object, permutations, ...)
{
    nperm <- nrow(permutations)
    ## Refuse to handle models with missing data
    if (!is.null(object$na.action))
        stop("by = 'margin' models cannot handle missing data")
    ## We need term labels but without Condition() terms
    trms <- drop.scope(object)
    trmlab <- trms[trms %in% attr(terms(object$terminfo),
                                      "term.labels")]
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
                     permutations, ...))
    ## Chande in df
    Df <- sapply(mods, function(x) x$df[2]) - dfbig
    ## F of change
    Chisq <- sapply(mods, function(x) x$chi[2]) - chibig
    Fstat <- (Chisq/Df)/(chibig/dfbig)
    ## Simulated F-values
    Fval <- sapply(mods, function(x) x$den)
    ## Had we an empty model we need to clone the denominator
    if (length(Fval) == 1)
        Fval <- matrix(Fval, nrow=nperm)
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

### Marginal test for axes

`anovacca.byaxis` <-
    function(object, permutations, model, parallel)
{
    nperm <- nrow(permutations)
    ## Observed F-values and Df
    eig <- object$CCA$eig
    resdf <- nobs(object) - length(eig) - max(object$pCCA$rank, 0) - 1
    Fstat <- eig/object$CA$tot.chi*resdf
    Df <- rep(1, length(eig))
    ## Marginal P-values
    LC <- object$CCA$u
    if (!is.null(object$na.action))
        LC <- napredict(structure(object$na.action,
                                  class="exclude"), LC)
    LC <- as.data.frame(LC)
    fla <- reformulate(names(LC))
    Pvals <- numeric(length(eig))
    environment(object$terms) <- environment()
    for (i in 1:length(eig)) {
        part <- paste("~ . +Condition(",
                      paste(names(LC)[-i], collapse = "+"), ")")
        upfla <- update(fla, part)
        ## only one axis, and cannot partial out?
        if (length(eig) == 1)
            mod <- permutest(object, permutations, model = model,
                             parallel = parallel)
        else
            mod <-
                permutest(update(object, upfla, data = LC),
                          permutations, model = model,
                          parellel = parallel)
        Pvals[i] <- (sum(mod$F.perm >= mod$F.0) + 1)/(nperm+1)
    }
    out <- data.frame(c(Df, resdf), c(eig, object$CA$tot.chi),
                      c(Fstat, NA), c(Pvals,NA))
    rownames(out) <- c(names(eig), "Residual")
    colnames(out) <- c("Df", "Chisq", "F", "Pr(>F)")
    class(out) <- c("anova", "data.frame")
    out
}
