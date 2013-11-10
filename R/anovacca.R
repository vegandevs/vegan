`anovacca` <-
    function(object, ..., permutations = how(nperm=999), by = NULL,
             strata = NULL) 
{
    if (is.null(object$CA) || is.null(object$CCA) ||
        object$CCA$rank == 0 || object$CA$rank == 0)
        return(anova.ccanull(object))
    if (!is.null(by)) {
        by <- match.arg(by, c("axis", "terms", "margin"))
        .NotYetUsed("by")
    }
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
        permutations <- shuffleSet(nrow(object$CCA$u),
                                   control = permutations)
    else # we got a permutation matrix and seed is unknown
        seed <- NA
    nperm <- nrow(permutations)
    ## stop permutations block
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
