`add1.cca`<-
    function(object, scope, test = c("none", "permutation"),
             permutations = how(nperm = 199), ...)
{
    if (inherits(object, "prc"))
        stop("'step'/'add1' cannot be used for 'prc' objects")
    if (is.null(object$terms))
        stop("ordination model must be fitted using formula")
    test <- match.arg(test)
    ## Default add1
    # don't show messages about aliased terms
    out <- suppressMessages(NextMethod("add1", object, test = "none"))
    cl <- class(out)
    ## Loop over terms in 'scope' and do anova.cca
    if (test == "permutation") {
        ## Avoid nested Condition(Condition(x) + z)
        if (!is.character(scope))
            scope <- add.scope(object, update.formula(object, scope))
        ns <- length(scope)
        adds <- matrix(0, ns+1, 2)
        adds[1, ] <- NA
        for (i in 1:ns) {
            tt <- scope[i]
            ## anova.ccalist for 2..scope terms
            if (!is.null(object$CCA)) {
                nfit <- suppressMessages(
                    update(object, as.formula(paste(". ~ . +", tt))))
                tmp <- anova(object, nfit, permutations = permutations, ...)
                adds[i+1,] <- unlist(tmp[2, 5:6])
            }
            else { # first term: simple anova
                nfit <- update(object,
                               as.formula(paste(". ~ . +", tt)))
                tmp <- anova(nfit,  permutations = permutations, ...)
                adds[i+1,] <- unlist(tmp[1, 3:4])
            }
        }
        out <- cbind(out, adds)
        colnames(out) <- c("Df", "AIC", "F", "Pr(>F)")
        ## check for redundant (0 Df) terms
        if (any(nas <- out[,1] < 1, na.rm = TRUE)) {
            out[[3]][nas] <- NA
            out[[4]][nas] <- NA
        }
        class(out) <- cl
    }
    out
}
