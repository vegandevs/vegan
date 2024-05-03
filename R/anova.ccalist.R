`anovaCCAlist` <-
    function(object, permutations, model, parallel)
{
    EPS <- sqrt(.Machine$double.eps)
    ## 'object' *must* be a list of cca objects, and 'permutations'
    ## *must* be a permutation matrix -- we assume that calling
    ## function takes care of this, and this function is not directly
    ## called by users.
    nmodels <- length(object)
    ## check that input is valid
    ## 1. All models must be fitted with the same method
    method <- sapply(object, function(z) z$method)
    if (!all(method == method[1]))
        stop("same ordination method must be used in all models")
    else
        method <- method[1]
    ## 2. All models must be fitted with formula interface
    if (any(sapply(object, function(x) is.null(x$terms))))
        stop("all models must be fitted with formula interface")
    ## 3. Same response
    resp <- sapply(object, function(z) deparse(formula(z)[[2]]))
    if (!all(resp == resp[1]))
        stop("response must be same in all models")
    ## 4. Same no. of observations
    N <- sapply(object, nobs)
    if (!all(N == N[1]))
        stop("number of observations must be same in all models")
    else
        N <- N[1]
    ## 5. Terms must be nested
    trms <- lapply(object, function(z) labels(terms(z)))
    o  <- order(sapply(trms, length))
    for (i in 2:nmodels)
        if (!all(trms[[o[i-1]]] %in% trms[[o[i]]]))
            stop("models must be nested")

    ## Check permutation matrix
    nperm <- nrow(permutations)
    ## check
    if (ncol(permutations) != N)
        stop(gettextf("'permutations' have %d columns, but data have %d rows",
                      ncol(nperm), N))
    ## All models are evaluated in permutest.cca with identical
    ## permutations so that the differences of single permutations can
    ## be used to assess the significance of differences of fitted
    ## models. This strictly requires nested models (not checked
    ## here): all terms of the smaller model must be included in the
    ## larger model.
    mods <- lapply(object, function(z)
                   permutest.cca(z, permutations = permutations,
                                 model = model, parallel = parallel))
    dfs <- sapply(mods, function(z) z$df)
    dev <- sapply(mods, function(z) z$chi)
    resdf <- dfs[2,]
    df <- -diff(resdf)
    resdev <- dev[2,]
    changedev <- -diff(resdev)
    big <- which.min(resdf)
    scale <- resdev[big]/resdf[big]
    fval <- changedev/df/scale
    ## Collect permutation results: denominator of F varies in each
    ## permutation.
    pscale <- mods[[big]]$den/resdf[big]
    ## Numerator of F
    pfvals <- sapply(mods, function(z) z$num)
    if (is.list(pfvals))
        pfvals <- do.call(cbind, pfvals)
    pfvals <- apply(pfvals, 1, diff)
    ## dropped to vector?
    if (!is.matrix(pfvals))
        pfvals <- matrix(pfvals, nrow = 1, ncol = nperm)
    pfvals <- sweep(pfvals, 1, df, "/")
    pfvals <- sweep(pfvals, 2, pscale, "/")
    pval <- rowSums(sweep(pfvals, 1, fval - EPS, ">="))
    pval <- (pval + 1)/(nperm + 1)
    ## collect table
    table <- data.frame(resdf, resdev, c(NA, df),
                        c(NA,changedev), c(NA,fval), c(NA,pval))
    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    dimnames(table) <- list(1L:nmodels,
                            c("ResDf", paste0("Res", varname), "Df",
                              varname, "F", "Pr(>F)"))
    ## Collect header information
    formulae <- sapply(object,
                       function(z) deparse(formula(z), width.cutoff = 500))
    head <- paste0("Permutation tests for ", method, " under ",
                  mods[[big]]$model, " model\n",
                   howHead(attr(permutations, "control")))
    topnote <- paste("Model ", format(1L:nmodels), ": ", formulae,
                     sep = "", collapse = "\n")
    structure(table, heading = c(head,topnote),
              F.perm = t(pfvals),
              class = c("anova.cca", "anova", "data.frame"))
}
