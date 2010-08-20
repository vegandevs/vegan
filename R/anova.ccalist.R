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
    resdf <- sapply(objects, function(x) x$CA$rank)
    resdev <- sapply(objects, function(x) x$CA$tot.chi)
    ## residual df (df) do not necessarily decrease with rank deficit
    ## models: we also need model (CCA) df and their change
    moddf <- sapply(objects, function(x) x$CCA$qrank)
    ## ANOVA table except test
    table <- data.frame(moddf, resdf, resdev, c(NA, diff(moddf)),
                        c(NA, -diff(resdev)))
    dimnames(table) <- list(1:nmodels, c("Model.Df", "Res.Df", "RSS", "Df", "Sum of Sq"))
    variables <- sapply(objects,
                        function(x) paste(deparse(formula(x)), collapse="\n"))
    title <- "Analysis of Variance Table\n"
    topnote <- paste("Model ", format(1:nmodels), ": ", variables,
                     sep = "", collapse = "\n")
    table <- structure(table, heading = c(title, topnote))
    ##Collect tests
    mods <- list()
    for(i in 1:nmodels) {
        if (i > 1)
            assign(".Random.seed", mods[[1]]$Random.seed,
                   envir = .GlobalEnv)
        mods[[i]] <- permutest(objects[[i]])
    }
    class(table) <- c("anova", class(table))
    table
}
