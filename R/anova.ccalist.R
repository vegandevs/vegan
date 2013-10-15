`anova.ccalist` <-
    function(object, ..., permutations = 99)
{
    ## Collect cca class objects. FIXME: Eventually this should be in
    ## a function that calls permutest.ccalist after collecting model
    ## objects from dotargs.
    dotargs <- list(...)
    if (length(dotargs)) {
        isCCA <- sapply(dotargs, function(z) inherits(z, "cca"))
        dotargs <- dotargs[isCCA]
        if (length(dotargs))
            object <- c(list(object), dotargs)
    }
    nmodels <- length(object)
    ## check that input is valid
    ## 1. All models must be fitted with the same method
    method <- sapply(object, function(z) z$method)
    if (!all(method == method[1]))
        stop("same ordination method must be used in all models")
    else
        method <- method[1]
    ## 2. Same response
    resp <- sapply(object, function(z) deparse(formula(z)[[2]]))
    if (!all(resp == resp[1]))
        stop("response must be same in all models")
    ## 3. Same no. of observations
    N <- sapply(object, nobs)
    if (!all(N = N[1]))
        stop("number of observations must be same in all models")
    else
        N <- N[1]
    ## 4. Terms must be nested
    trms <- lapply(object, function(z) labels(terms(z)))
    o  <- order(sapply(trms, length))
    for(i in 2:nmodels) 
        if(!all(trms[[o[i-1]]] %in% trms[[o[i]]]))
            stop("models must be nested")
        
    ## Create permutation matrix if it does not exist. FIXME: should
    ## take arguments for restricted permutation
    if (length(permutations) == 1) 
            permutations <- shuffleSet(N, permutations)
    ## permutations is now matrix
    nperm <- nrow(permutations)
    ## check
    if (ncol(permutations) != N)
        stop(gettextf("permutation matrix has %d columns, but you have %d sites",
                      ncol(nperm), N))
    ## All models are evaluated in permutest.cca with identical
    ## permutations so that the differences of single permutations can
    ## be used to assess the significance of differences of fitted
    ## models. This strictly requires nested models (not checked
    ## here): all terms of the smaller model must be included in the
    ## larger model. FIXME: should pass arguments to permutest.cca.
    mods <- lapply(object, function(z)
                   permutest(z, permutations = permutations))
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
        pfvals <- matrix(pfvals, nrow=1, ncol=nperm)
    pval <- rowSums(sweep(pfvals, 1, fval, ">="))
    pval <- (pval + 1)/(nperm+1)
    pfvals <- sweep(pfvals, 1, df, "/")
    pfvals <- sweep(pfvals, 2, pscale, "/")
    ## collect table
    table <- data.frame(resdf, resdev, c(NA, df),
                        c(NA,changedev), c(NA,fval), c(NA,pval))
    dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Res. Chisq", 
                                          "Df", "Chisq", "F", "Pr(>F)"))
    ## Collect header information
    formulae <- sapply(object, function(z) deparse(formula(z)))
    head <- paste0("Permutation tests for ", method, " under ",
                  mods[[big]]$model, " model\nwith ", nperm,
                   " permutations\n")
    topnote <- paste("Model ", format(1L:nmodels), ": ", formulae,
                     sep = "", collapse = "\n")
    structure(table, heading=c(head,topnote), class = c("anova", "data.frame"))
}
