### anova.cca cannot be performed if residuals or constraints are
### NULL, and this function handles these cases (but it doesn't test
### that these are the cases).
`anova.ccanull` <-
    function(object, ...)
{
    table <- matrix(0, nrow = 2, ncol = 4)
    if (object$CA$rank == 0) {
        table[1,] <- c(object$CCA$qrank, object$CCA$tot.chi, NA, NA)
        table[2,] <- c(0,0,NA,NA)
    }
    else {
        table[1,] <- c(0,0,0,NA)
        table[2,] <- c(nrow(object$CA$u) - 1, object$CA$tot.chi, NA, NA)
    }
    rownames(table) <- c("Model", "Residual")
    if (inherits(object, c("capscale", "dbrda")) &&
        (object$adjust != 1 || is.null(object$adjust)))
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    colnames(table) <- c("Df", varname, "F", "Pr(>F)")
    table <- as.data.frame(table)
    if (object$CA$rank == 0)
        head <- "No residual component\n"
    else if (is.null(object$CCA) || object$CCA$rank == 0)
        head <- "No constrained component\n"
    else
        head <- c("!!!!! ERROR !!!!!\n")
    head <- c(head, paste("Model:", c(object$call)))
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    else
        seed <- NULL
    structure(table, heading = head, Random.seed = seed,
              class = c("anova.cca", "anova", "data.frame"))
}
