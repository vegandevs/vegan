### Forward selection to maximize R2.adjusted, but stopping once the
### R2.adjusted of the maximum model ('scope') is exceeded, after
### Blanchet, Legendre & Borcard: Ecology 89, 2623--2623; 2008.

`ordiR2step` <-
    function(object, scope, direction = c("both", "forward"),
             Pin = 0.05, R2scope = TRUE, permutations = how(nperm=499),
             trace = TRUE, ...)
{
    direction <- match.arg(direction)
    if (missing(scope))
        stop("needs scope")
    ## Works only for rda(): cca() does not have (yet) R2.adjusted
    if (!inherits(object, "rda"))
        stop("can be used only with rda() or capscale()")
    ## No R2 for capscale with negative eigenvalues
    if (inherits(object, "capscale") && !is.null(object$CA$imaginary.chi))
        stop("cannot be used when capscale() has negative eigenvalues")
    ## Get R2 of the original object
    if (is.null(object$CCA))
        R2.0 <- 0
    else
        R2.0 <- RsquareAdj(object)$adj.r.squared
    ## Get R2 of the scope
    if (inherits(scope, "rda")) 
        scope <- delete.response(formula(scope))
    if (!inherits(scope, "formula"))
        scope <- reformulate(scope)
    if (R2scope)
        R2.all <- RsquareAdj(update(object, scope))
    else
        R2.all <- list(adj.r.squared = NA)
    ## Check that the full model can be evaluated
    if (is.na(R2.all$adj.r.squared) && R2scope) 
        stop("the upper scope cannot be fitted (too many terms?)")
    R2.all <- R2.all$adj.r.squared
    ## Collect data to anotab returned as the 'anova' object
    anotab <-  list()
    ## Step forward and continue as long as R2.adj improves and R2.adj
    ## remains below R2.adj < R2.all
    R2.previous <- R2.0
    drops <- NULL
    repeat {
        if (trace) {
            cat("Step: R2.adj=", R2.previous, "\n")
            cat(pasteCall(formula(object)), "\n")
        }
        adds <- add.scope(object, scope)
        if (direction == "both")
            drops <- drop.scope(object)
        ## Nothing to add or drop, and we're done: break
        if (length(adds) == 0 && length(drops) == 0)
            break
        R2.adds <- numeric(length(adds) + length(drops))
        if (length(adds))
            adds <- paste("+", adds)
        if (length(drops))
            drops <- paste("-", drops)
        names(R2.adds) <- c(adds, drops) 
        ## Loop over add scope
        for (trm in seq_along(R2.adds)) {
            fla <- paste(". ~ .", names(R2.adds[trm]))
            R2.tmp <- RsquareAdj(update(object, fla))$adj.r.squared
            if (!length(R2.tmp))
                R2.tmp <- 0
            R2.adds[trm] <- R2.tmp
        }
        best <- which.max(R2.adds)
        if (trace) {
            out <- sort(c("<All variables>" = R2.all, "<none>" = R2.previous,
                          R2.adds), decreasing = TRUE)
            out <- as.matrix(out)
            colnames(out) <- "R2.adjusted"
            print(out)
            cat("\n")
        }
        ## See if the best should be kept
        ## First criterion: R2.adj improves and is still lower or
        ## equal than for the full model of the scope
        if (R2.adds[best] > R2.previous &&
            if (R2scope) R2.adds[best] <= R2.all else TRUE) {
            ## Second criterion: added variable is significant
            tst <- add1(object, scope = adds[best], test="permu",
                        permutations = permutations,
                        alpha = Pin, trace = FALSE, ...)
            if (trace) {
                print(tst[-1,])
                cat("\n")
            }
            if (tst[,"Pr(>F)"][2] > Pin)
                break
            fla <- paste("~  .", adds[best])
            object <- update(object, fla)
            R2.previous <- RsquareAdj(object)$adj.r.squared
            anotab <- rbind(anotab, cbind("R2.adj" = R2.previous, tst[2,]))
        } else {
            break
        }
    }
    if (NROW(anotab) > 0) {
        anotab <- rbind(anotab, "<All variables>" = c(R2.all, rep(NA, 4)))
        class(anotab) <- c("anova", class(anotab))
        object$anova <- anotab
    }
    object
}
