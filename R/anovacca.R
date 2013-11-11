`anovacca` <-
    function(object, ..., permutations = how(nperm=999), by = NULL,
             strata = NULL) 
{
    if (!exists(".Random.seed", envir = .GlobalEnv,
                inherits = FALSE)) 
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    ## permutations is either a single number, a how() structure or a
    ## permutation matrix
    if (length(permutations) == 1) {
        nperm <- permutations
        permutations <- how(nperm = nperm)
    }
    if (!is.null(strata)) {
        if (!inherits(permutations, "how"))
            stop("'strata' can be used only with simple permutation or with 'how()'")
        if (!is.null(permutations$block))
            stop("'strata' cannot be applied when 'blocks' are defined in 'how()'")
        permutations <- update(permutations, blocks = strata)
    }
    ## now permutations is either a how() structure or a permutation
    ## matrix. Make it to a matrix if it is "how"
    if (inherits(permutations, "how"))
        permutations <- shuffleSet(nrow(object$CA$u),
                                   control = permutations)
    else # we got a permutation matrix and seed is unknown
        seed <- NA
    nperm <- nrow(permutations)
    ## stop permutations block
    ## see if this was a list of ordination objects
    dotargs <- list(...)
    if (length(dotargs)) {
        isCCA <- sapply(dotargs, function(z) inherits(z, "cca"))
        if (any(isCCA)) {
            ## we do not want to give dotargs to anova.ccalist, but we
            ## evaluate 'parallel' and 'model' here
            if (is.null(dotargs$model))
                model <- "reduced"
            else
                model <- dotargs$model
            if (is.null(dotargs$parallel))
                parallel <- NULL
            else
                parallel <- dotargs$parallel
            dotargs <- dotargs[isCCA]
            object <- c(list(object), dotargs)
            sol <-
                anova.ccalist(object, 
                              permutations = permutations,
                              model = model,
                              parallel = parallel)
            return(sol)
        }
    }
    ## by cases
    if (!is.null(by)) {
        by <- match.arg(by, c("axis", "terms", "margin"))
        .NotYetUsed("by")
    }
    ## basic overall test
    tst <- permutest.cca(object, permutations = permutations, ...)
    Fval <- c(tst$F.0, NA)
    Pval <- (sum(tst$F.perm >= tst$F.0) + 1)/(tst$nperm + 1)
    Pval <- c(Pval, NA)
    nperm <- c(tst$nperm, NA)
    table <- data.frame(tst$df, tst$chi, Fval, nperm, Pval)
    is.rda <- inherits(object, "rda")
    colnames(table) <- c("Df", ifelse(is.rda, "Var", "Chisq"), 
                         "F", "N.Perm", "Pr(>F)")
    head <- paste("Permutation test for", tst$method, "under", 
                  tst$model, "model\n")
    if (!is.null(tst$strata)) 
        head <- paste(head, "Permutations stratified within '", 
                      tst$strata, "'\n", sep = "")
    mod <- paste("Model:", c(object$call))
    structure(table, heading = c(head, mod), Random.seed = seed, 
              class = c("anova.cca", "anova", "data.frame"))
}
