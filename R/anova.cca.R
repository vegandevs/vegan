`anova.cca` <-
    function(object, ..., permutations = how(nperm=999), by = NULL,
             model = c("reduced", "direct", "full"),
             parallel = getOption("mc.cores"), strata = NULL,
             cutoff = 1, scope = NULL, test = c("permutation", "F"),
             statistic = c("pseudo-F", "modified-F"), groups = NULL)
{
    EPS <- sqrt(.Machine$double.eps) # for permutation P-values
    model <- match.arg(model)
    test <- match.arg(test)
    statistic <- match.arg(statistic)
    if (test == "F") {
        if (statistic == "modified-F")
            stop("'statistic = \"modified-F\"' requires a permutation test")
        permutations <- 0
    }
    ## permutation matrix
    N <- nrow(object$CA$u)
    permutations <- getPermuteMatrix(permutations, N, strata = strata)
    seed <- attr(permutations, "seed")
    control <- attr(permutations, "control")
    ## see if this was a list of ordination objects
    dotargs <- list(...)
    ## we do not want to give dotargs to anovaCCAlist, but we
    ## evaluate 'parallel' and 'model' here
    if (length(dotargs)) {
        isCCA <- sapply(dotargs, function(z) inherits(z, "cca"))
        if (any(isCCA)) {
            if (statistic == "modified-F")
                stop("'statistic = \"modified-F\"' cannot compare multiple models")
            dotargs <- dotargs[isCCA]
            object <- c(list(object), dotargs)
            sol <- anovaCCAlist(object,
                                permutations = permutations,
                                model = model,
                                parallel = parallel)
            attr(sol, "Random.seed") <- seed
            attr(sol, "control") <- control
            if (test == "F") {
                big <- which.max(sapply(object, function(z)
                    length(labels(terms(z)))))
                sol <- anovaCCAparametric(object[[big]], sol)
            }
            return(sol)
        }
    }
    ## We only have a single model: check if it is empty
    if (is.null(object$CA) || is.null(object$CCA) ||
        object$CCA$rank == 0 || object$CA$rank == 0) {
        if (statistic == "modified-F")
            stop("the modified-F test requires non-empty constrained and ",
                 "residual spaces")
        return(anovaCCAnull(object))
    }
    ## by cases
    if (!is.null(by)) {
        by <- match.arg(by, c("terms", "margin", "axis", "onedf"))
        if (statistic == "modified-F")
            stop("'by' must be NULL with 'statistic = \"modified-F\"'")
        if (by %in% c("terms", "margin") && is.null(object$terms))
            stop("model must be fitted with formula interface")
        sol <- switch(by,
                      "terms" = anovaCCAbyterm(object,
                      permutations = permutations,
                      model = model, parallel = parallel),
                      "margin" = anovaCCAbymargin(object,
                      permutations = permutations,
                      model = model, parallel = parallel,
                      scope = scope),
                      "axis" = anovaCCAbyaxis(object,
                      permutations = permutations,
                      model = model, parallel = parallel,
                      cutoff = cutoff),
                      "onedf" = anovaCCAby1df(object,
                       permutations = permutations,
                       model = model, parallel = parallel)
                      )
        attr(sol, "Random.seed") <- seed
        attr(sol, "control") <- control
        if (test == "F")
            sol <- anovaCCAparametric(object, sol)
        return(sol)
    }
    ## basic overall test: pass other arguments except 'strata'
    ## because 'permutations' already is a permutationMatrix
    tst <- permutest.cca(object, permutations = permutations,
                         model = model, parallel = parallel, ...,
                         statistic = statistic, groups = groups)
    Fval <- c(tst$F.0, NA)
    Pval <- (sum(tst$F.perm >= tst$F.0 - EPS) + 1)/(tst$nperm + 1)
    Pval <- c(Pval, NA)
    table <- data.frame(tst$df, tst$chi, Fval, Pval)
    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    colnames(table) <- c("Df", varname, "F", "Pr(>F)")
    stathead <- if (statistic == "modified-F")
        "Modified pseudo-F allowing heterogeneous dispersions\n"
    else ""
    head <- paste0("Permutation test for ", tst$method, " under ",
                  tst$model, " model\n", stathead, howHead(control))
    mod <- paste("Model:", c(object$call))
    if (varname == "SumOfSqs")
        mod <- paste0(mod, "\nSumOfSqs is ", object$inertia)
    sol <- structure(table, heading = c(head, mod), Random.seed = seed,
                     control = control, F.perm = tst$F.perm,
                     class = c("anova.cca", "anova", "data.frame"))
    if (statistic == "modified-F") {
        attr(sol, "statistic") <- statistic
        attr(sol, "denominator") <- tst$den.0
        attr(sol, "denominator.df") <- tst$den.df
        attr(sol, "residual.df") <- tst$res.df
        attr(sol, "groups") <- tst$groups
        attr(sol, "group.size") <- tst$group.size
        attr(sol, "group.dispersion") <- tst$group.dispersion
    }
    if (test == "F")
        sol <- anovaCCAparametric(object, sol)
    sol
}
