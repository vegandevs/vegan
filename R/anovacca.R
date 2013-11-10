`anovacca` <-
    function(object, ..., permutations = how(nperm=999), by = NULL) 
{
    if (is.null(object$CA) || is.null(object$CCA) ||
        object$CCA$rank == 0 || object$CA$rank == 0)
        return(anova.ccanull(object))
    if (!is.null(by)) {
        by <- match.arg(by, c("axis", "terms", "margin"))
        .NotYetUsed("by")
    }
    seed <- NULL
    tst <- permutest.cca(object, permutations = permutations, ...)
    if (is.null(seed)) 
        seed <- tst$Random.seed
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
