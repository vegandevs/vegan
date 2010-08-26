### Forward selection to maximize R2.adjusted, but stopping once the
### R2.adjusted of the maximum model ('scope') is exceeded, after
### Blanchet, Legendre & Borcard: Ecology 89, 2623--2623; 2008.

`ordiR2step` <-
    function(object, scope, Pin = 0.05, pstep = 100,
             perm.max = 1000, trace = TRUE, ...)
{
    if (missing(scope))
        stop("needs scope")
    ## Works only for rda(): cca() does not have (yet) R2.adjusted
    if (!inherits(object, "rda"))
        stop("can be used only with rda() or capscale()")
    ## Get R2 of the original object
    if (is.null(object$CCA))
        R2.0 <- 0
    else
        R2.0 <- RsquareAdj(object)$adj.r.squared
    ## Get R2 of the scope
    if (inherits(scope, "rda")) {
        R2.all <- RsquareAdj(scope)$adj.r.squared
        scope <- formula(scope)
    } else {
        if (!inherits(scope, "formula"))
            scope <- reformulate(scope)
        R2.all <- RsquareAdj(update(object, scope))$adj.r.squared
    }
    ## Check that the full model can be evaluated
    if (is.na(R2.all))
        stop("the upper scope cannot be fitted (too many terms?)")
    ## Step forward and continue as long as R2.adj improves and R2.adj
    ## remains below R2.adj < R2.all
    R2.previous <- R2.0
    repeat {
        if (trace) {
            cat("Step: R2.adj=", R2.previous, "\n")
            cat(pasteCall(formula(object)), "\n")
        }
        adds <- add.scope(object, scope)
        ## Nothing to add, and we're done: break
        if (length(adds) == 0)
            break
        R2.adds <- numeric(length(adds))
        names(R2.adds) <- adds
        ## Loop over add scope
        for (trm in seq_along(adds)) {
            fla <- paste("~  . + ", adds[trm])
            R2.adds[trm] <- RsquareAdj(update(object, fla))$adj.r.squared
        }
        best <- which.max(R2.adds)
        if (trace) {
            names(R2.adds) <- paste("+", names(R2.adds))
            out <- sort(c("<All variables>" = R2.all, R2.adds), decreasing = TRUE)
            out <- as.matrix(out)
            colnames(out) <- "R2.adjusted"
            print(out)
            cat("\n")
        }
        ## See if the best should be kept
        ## First criterion: R2.adj improves and is still lower or
        ## equal than for the full model of the scope
        if (R2.adds[best] > R2.previous && R2.adds[best] <= R2.all) {
            ## Second criterion: added variable is significant
            tst <- add1(object, scope = adds[best], test="permu",
                        pstep = pstep, perm.max = perm.max,
                        alpha = Pin, trace = FALSE, ...)
            if (trace) {
                print(tst[-1,])
                cat("\n")
            }
            if (tst[,"Pr(>F)"][2] > Pin)
                break
            fla <- paste("~  . +", adds[best])
            object <- update(object, fla)
            R2.previous <- RsquareAdj(object)$adj.r.squared
        } else {
            break
        }
    }
    object
}
