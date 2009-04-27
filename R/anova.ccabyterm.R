`anova.ccabyterm` <-
    function (object, step = 100, ...) 
{
    trm <- terms(object)
    call <- paste("Model:", c(object$call))
    trmlab <- attr(trm, "term.labels")
    trmlab <- trmlab[trmlab %in% attr(terms(object$terminfo), 
                                      "term.labels")]
    ntrm <- length(trmlab)
    ## 'adj' puts the result together with the permutations and reduces
    ## number of simulations by one so that P = (hits+1)/(permutations+1).
    ## The first step is reduced by adj.
    adj <- (step %% 10) == 0
    step <- step - adj
    pchi <- matrix(0, nrow = ntrm + 1, ncol = step)
    chi <- numeric(ntrm + 1)
    df <- numeric(ntrm + 1)
    names(df) <- c(trmlab, "Residual")
    sim <- permutest.cca(object, permutations = step, ...)
    pchi[ntrm + 1, ] <- sim$den
    pchi[ntrm, ] <- sim$num
    df[ntrm:(ntrm + 1)] <- sim$df
    chi[ntrm:(ntrm + 1)] <- sim$chi
    for (.ITRM in ntrm:2) {
        if (ntrm < 2) 
            break
        assign(".Random.seed", sim$Random.seed, envir = .GlobalEnv)
        fla <- as.formula(paste(" . ~ . -", trmlab[.ITRM]))
        object <- update(object, fla)
        if (is.null(object$CCA)) 
            break
        sim <- permutest.cca(object, permutations = step, ...)
        pchi[.ITRM, ] <- pchi[.ITRM, ] - sim$num
        chi[.ITRM] <- chi[.ITRM] - sim$chi[1]
        df[.ITRM] <- df[.ITRM] - sim$df[1]
        pchi[.ITRM - 1, ] <- sim$num
        chi[.ITRM - 1] <- sim$chi[1]
        df[.ITRM - 1] <- sim$df[1]
    }
    Fval <- chi/df/(chi[ntrm + 1]/df[ntrm + 1])
    Fval[ntrm + 1] <- NA
    pchi <- sweep(pchi, 1, df, "/")
    pchi[-(ntrm + 1), ] <- sweep(pchi[-(ntrm + 1), , drop = FALSE], 
                                 2, pchi[ntrm + 1, , drop = FALSE], "/")
    P <- rowSums(sweep(pchi[-(ntrm + 1), , drop = FALSE], 1, 
                       Fval[-(ntrm + 1)], ">"))
    P <- c((P + adj)/(step + adj), NA)
    out <- data.frame(df, chi, Fval, c(rep(step, ntrm), NA), 
                      P)
    inertname <- if (sim$method == "cca") 
        "Chisq"
    else "Var"
    colnames(out) <- c("Df", inertname, "F", "N.Perm", "Pr(>F)")
    out <- out[out[, 1] > 0 | out[, 2] > sqrt(.Machine$double.eps), 
               ]
    head <- paste("Permutation test for", sim$method, "under", 
                  sim$model, "model\nTerms added sequentially (first to last)\n")
    if (!is.null(sim$strata)) 
        head <- paste(head, "Permutations stratified within `", 
                      sim$strata, "'\n", sep = "")
    structure(out, heading = c(head, call), Random.seed = sim$Random.seed, 
              class = c("anova.cca", "anova", "data.frame"))
}
