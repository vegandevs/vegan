`anova.cca` <-
    function (object, alpha = 0.05, beta = 0.01, step = 100, perm.max = 9999, 
              by = NULL, ...) 
{
    if (is.null(object$CA) || is.null(object$CCA))
        return(anova.ccanull(object))
    perm.max <- max(step-1, perm.max)
    if (perm.max %% step == 0)
        perm.max <- perm.max - 1
    if (!is.null(by)) {
        by <- match.arg(by, c("axis", "terms", "margin"))
        if (by == "axis") 
            sol <- anova.ccabyaxis(object, alpha = alpha, beta = beta, 
                                   step = step, perm.max = perm.max, by = NULL, 
                                   ...)
        else if (by == "margin") {
            sol <- anova.ccabymargin(object, alpha = alpha, beta = beta,
                                     step = step, perm.max = perm.max,
                                     by = NULL, ...)
            }
        else {
            mf <- match.call(expand.dots = FALSE)
            if (!is.null(mf$...) && any(k <- pmatch(names(mf$...), 
                                                    "permutations", nomatch = FALSE))) 
                step <- unlist(mf$...[k == 1])
            sol <- anova.ccabyterm(object, step = step, ...)
        }
        return(sol)
    }
    seed <- NULL
    betaq <- c(beta/2, 1 - beta/2)
    nperm <- 0
    unsure <- TRUE
    hits <- 0
    while (unsure && nperm < perm.max) {
        adj <- as.numeric(nperm == 0)
        tst <- permutest.cca(object, step - adj, ...)
        if (is.null(seed)) 
            seed <- tst$Random.seed
        nperm <- nperm + step - adj
        hits <- hits + sum(tst$F.perm >= tst$F.0)
        fork <- qbinom(betaq, nperm, alpha)
        if (hits < fork[1] || hits > fork[2]) 
            unsure <- FALSE
    }
    Fval <- c(tst$F.0, NA)
    Pval <- c((hits+1)/(nperm+1), NA)
    nperm <- c(nperm, NA)
    table <- data.frame(tst$df, tst$chi, Fval, nperm, Pval)
    is.rda <- inherits(object, "rda")
    colnames(table) <- c("Df", ifelse(is.rda, "Var", "Chisq"), 
                         "F", "N.Perm", "Pr(>F)")
    head <- paste("Permutation test for", tst$method, "under", 
                  tst$model, "model\n")
    if (!is.null(tst$strata)) 
        head <- paste(head, "Permutations stratified within `", 
                      tst$strata, "'\n", sep = "")
    mod <- paste("Model:", c(object$call))
    structure(table, heading = c(head, mod), Random.seed = seed, 
              class = c("anova.cca", "anova", "data.frame"))
}
