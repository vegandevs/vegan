### Parametric anova.cca: using F-statistic with Greenhouse-Geisser
### asphericity correction

### input: ordination 'object' and anova.cca table 'anotab'. Replaces
### P-values with the parametric estimate from F-distribution.

`anovaCCAparametric` <-
    function(object, anotab, ...)
{
    if (!inherits(object, "rda") || inherits(object, c("dbrda", "capscale")))
        stop("parametric test only implemented for rda")
    ## Greenhouse-Geisser epsilon
    NC <- nrow(object$CCA$v)
    ssd <- SSD(object, type = "response")$SSD
    u <- solve(diag(nrow = NC), ssd)
    lam <- Re(eigen(u, only.values = TRUE)$values)
    GG.eps <- sum(lam)^2 / sum(lam^2) / NC
    ## Parametric F-values. anovaCCAlist and other methods have
    ## different numbers of columns and residual Df in different place.
    if (ncol(anotab) == 4) {
        dfres <- anotab[nrow(anotab), "Df"] * NC * GG.eps
        df <- anotab[-nrow(anotab), "Df"] * NC * GG.eps
    } else if (ncol(anotab) == 6) { # anovaCCAlist
        dfres <- anotab[nrow(anotab), "ResDf"] * NC * GG.eps
        df <- anotab[, "Df"] * NC * GG.eps
    }
    anotab[,"Pr(>F)"] <- pf(anotab[, "F"], df, dfres, lower.tail = FALSE)
    ## Edit heading
    head <- attr(anotab, "heading")[1]
    headlines <- strsplit(head, "\n", fixed = TRUE)[[1]]
    if (length(headlines) == 4)
        bycase <- paste0("\n", headlines[2])
    else
        bycase <- NULL
    head <- paste("Parametric test for", object$method,
                  bycase,
                  "\nGreenhouse-Geisser epsilon:",
                  format(GG.eps, digits = 4L), "\n")
    attr(anotab, "heading")[1] <- head
    anotab
}
