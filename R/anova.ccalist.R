### A test of concept function for significance test of two (or more)
### cca, rda or capscale models. A "test of concept" means that there
### are no sanity chekcs, but only the basic code with assumption that
### input is correct.

### The idea is to run permutest with the same random number seed and
### take the difference of permutation results (denominator and
### numerator of pseudo-F) and take P-values of that.  Assume the two
### models are m1 <- cca(y ~ x1) and m2 <- cca(y ~ x1 + x2). Currently
### the difference can be analysed as a partial model, where the
### common part is partialled out m12 <- cca(y ~ Condition(x1) +
### x2). The direct difference of two models suggested here does not
### produce the same results except with model = "direct", because the
### partial model will permute residuals after conditions.
`anova.ccalist` <- function(object, ...)
{
    objects <- list(object, ...)
    nmodels <- length(objects)
    ## Collect statistics
    N <- nrow(object$CA$u)
    ranks <- sapply(objects, function(x)
                    if (is.null(x$CCA$qrank)) 0 else x$CCA$qrank)
    resdf <- N - ranks - 1
    resdev <- sapply(objects, deviance)
    moddev <- c(NA, -diff(resdev))
    moddf <- c(NA, diff(ranks))
    ##Collect tests
    mods <- list()
    for(i in 1:nmodels) {
        if (i > 1)
            assign(".Random.seed", mods[[1]]$Random.seed,
                   envir = .GlobalEnv)
        mods[[i]] <- permutest(objects[[i]])
    }
    ## Differences of permutations. In permutation F values, numerator
    ## is taken from each model, but all use the same denominator from
    ## the largest model.
    bigmodel <- which.min(resdf)
    F <- moddev/moddf/resdev[bigmodel]*resdf[bigmodel]
    den <- mods[[bigmodel]]$den/mods[[bigmodel]]$df[2]
    Pval <- rep(NA, nmodels)
    for (i in 2:nmodels) {
        F.perm <- abs(mods[[i]]$num - mods[[i-1]]$num)/abs(moddf[i])/abs(den)
        Pval[i] <- (sum(F.perm >= F[i]) + 1)/(length(F.perm) + 1)
    }   
    ## ANOVA table 
    table <- data.frame(resdf, resdev, moddf, moddev, F, Pval)
    dimnames(table) <- list(1:nmodels, c("Res.Df", "RSS", "Df", "Sum of Sq",
                                         "F", "Pr(>F)"))
    variables <- sapply(objects,
                        function(x) paste(deparse(formula(x)), collapse="\n"))
    title <- "Analysis of Variance Table\n"
    topnote <- paste("Model ", format(1:nmodels), ": ", variables,
                     sep = "", collapse = "\n")
    table <- structure(table, heading = c(title, topnote))
    class(table) <- c("anova.cca", "anova", class(table))
    table
}
