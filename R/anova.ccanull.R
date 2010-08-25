### anova.cca cannot be performed if residuals or constraints are
### NULL, and this function handles these cases (but it doesn't test
### that these are the cases).
`anova.ccanull` <-
    function(object, ...)
{
    table <- matrix(0, nrow = 2, ncol = 5)
    if (is.null(object$CA)) {
        table[1,] <- c(object$CCA$qrank, object$CCA$tot.chi, NA, 0, NA)
        table[2,] <- c(0,0,NA,NA,NA)
    }
    else {
        table[1,] <- c(0,0,0,0,NA)
        table[2,] <- c(nrow(object$CA$u) - 1, object$CA$tot.chi, NA, NA, NA)
    }
    rownames(table) <- c("Model", "Residual")
    colnames(table) <-  c("Df",
                          if (inherits(object, "rda")) "Var" else "Chisq", 
                          "F", "N.Perm", "Pr(>F)")
    table <- as.data.frame(table)
    if (is.null(object$CA))
        head <- "No residual component\n"
    else if (is.null(object$CCA))
        head <- "No constrained component\n"
    else
        head <- c("!!!!! ERROR !!!!!\n")
    head <- c(head, paste("Model:", c(object$call)))
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    structure(table, heading = head, Random.seed = seed,
              class = c("anova.cca", "anova", "data.frame"))
}
