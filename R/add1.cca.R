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
    out <- suppressMessages(NextMethod("add1", object, test = "none", ...))
    cl <- class(out)
    ## Loop over terms in 'scope' and do anova.cca
    if (test == "permutation") {
        ## Avoid nested Condition(Condition(x) + z)
        hasfla <- update(terms(object$terminfo), . ~  Condition(.))
        if (!is.character(scope)) 
            scope <- add.scope(object, update.formula(object, scope))
        ns <- length(scope)
        adds <- matrix(0, ns+1, 2)
        adds[1, ] <- NA
        for (i in 1:ns) {
            tt <- scope[i]
            ## Condition(.) previous terms (if present)
            if (!is.null(object$CCA)) {
                fla <- update(hasfla, paste("~ . +", tt))
                nfit <- update(object, fla)
            }
            else
                nfit <- update(object,
                               as.formula(paste(". ~ . +", tt)))
            tmp <- anova(nfit,  permutations = permutations, ...)
            adds[i+1,] <- unlist(tmp[1,3:4])
        }
        colnames(adds) <- colnames(tmp)[3:4]
        out <- cbind(out, adds)
        ## check for redundant (0 Df) terms
        if (any(nas <- out[,1] < 1, na.rm = TRUE)) {
            out[[3]][nas] <- NA
            out[[4]][nas] <- NA
        }
        class(out) <- cl
    }
    out
}
