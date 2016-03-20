`anova.cca` <-
    function(object, ..., permutations = how(nperm=999), by = NULL,
             model = c("reduced", "direct", "full"),
             parallel = getOption("mc.cores"), strata = NULL,
             cutoff = 1, scope = NULL)
{
    EPS <- sqrt(.Machine$double.eps) # for permutation P-values
    model <- match.arg(model)
    ## permutation matrix
    N <- nrow(object$CA$u)
    permutations <- getPermuteMatrix(permutations, N, strata = strata)
    seed <- attr(permutations, "seed")
    control <- attr(permutations, "control")
    nperm <- nrow(permutations)
    ## see if this was a list of ordination objects
    dotargs <- list(...)
    ## we do not want to give dotargs to anova.ccalist, but we
    ## evaluate 'parallel' and 'model' here
    if (length(dotargs)) {
        isCCA <- sapply(dotargs, function(z) inherits(z, "cca"))
        if (any(isCCA)) {
            dotargs <- dotargs[isCCA]
            object <- c(list(object), dotargs)
            sol <-
                anova.ccalist(object,
                              permutations = permutations,
                              model = model,
                              parallel = parallel)
            attr(sol, "Random.seed") <- seed
            attr(sol, "control") <- control
            return(sol)
        }
    }
    ## We only have a single model: check if it is empty
    if (is.null(object$CA) || is.null(object$CCA) ||
        object$CCA$rank == 0 || object$CA$rank == 0)
        return(anova.ccanull(object))
    ## by cases
    if (!is.null(by)) {
        by <- match.arg(by, c("terms", "margin", "axis"))
        if (is.null(object$terms))
            stop("model must be fitted with formula interface")
        sol <- switch(by,
                      "terms" = anova.ccabyterm(object,
                      permutations = permutations,
                      model = model, parallel = parallel),
                      "margin" = anova.ccabymargin(object,
                      permutations = permutations,
                      model = model, parallel = parallel,
                      scope = scope),
                      "axis" = anova.ccabyaxis(object,
                      permutations = permutations,
                      model = model, parallel = parallel,
                      cutoff = cutoff))
        attr(sol, "Random.seed") <- seed
        attr(sol, "control") <- control
        return(sol)
    }
    ## basic overall test: pass other arguments except 'strata'
    ## because 'permutations' already is a permutationMatrix
    tst <- permutest.cca(object, permutations = permutations,
                         model = model, parallel = parallel, ...)
    Fval <- c(tst$F.0, NA)
    Pval <- (sum(tst$F.perm >= tst$F.0 - EPS) + 1)/(tst$nperm + 1)
    Pval <- c(Pval, NA)
    table <- data.frame(tst$df, tst$chi, Fval, Pval)
    if (inherits(object, c("capscale", "dbrda")) &&
        (object$adjust != 1 || is.null(object$adjust)))
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    colnames(table) <- c("Df", varname, "F", "Pr(>F)")
    head <- paste0("Permutation test for ", tst$method, " under ",
                  tst$model, " model\n", howHead(control))
    mod <- paste("Model:", c(object$call))
    structure(table, heading = c(head, mod), Random.seed = seed,
              control = control, F.perm = tst$F.perm,
              class = c("anova.cca", "anova", "data.frame"))
}
