### Forward selection to maximize R2.adjusted, but stopping once the
### R2.adjusted of the maximum model ('scope') is exceeded, after
### Blanchet, Legendre & Borcard: Ecology 89, 2623--2623; 2008.

`ordiR2step` <-
    function(object, scope, Pin = 0.05, R2scope = TRUE,
             permutations = how(nperm=499),
             trace = TRUE, R2permutations = 1000, ...)
{
    if (is.null(object$terms))
        stop("ordination model must be fitted using formula")
    if (missing(scope))
        stop("needs scope")
    if (inherits(scope, "cca"))
        scope <- delete.response(formula(scope))
    if (!inherits(scope, "formula"))
        scope <- reformulate(scope)
    ## Get R2 of the original object
    if (is.null(object$CCA))
        R2.0 <- 0
    else
        R2.0 <- RsquareAdj(object,
                           permutations = R2permutations, ...)$adj.r.squared
    ## only accepts upper scope
    if (is.list(scope) && length(scope) <= 2L)
        scope <- scope$upper
    if (is.null(scope) || !length(add.scope(object, scope)))
        stop("needs upper 'scope': no terms can be added")
    ## Get R2 of the scope
    if (R2scope)
        R2.all <- suppressMessages(
            RsquareAdj(update(object, delete.response(formula(scope))),
                             permutations = R2permutations, ...))
    else
        R2.all <- list(adj.r.squared = NA)
    ## Check that the full model can be evaluated
    if (is.na(R2.all$adj.r.squared) && R2scope)
        stop("R2 cannot be adjusted in upper scope (too many terms?)")
    R2.all <- R2.all$adj.r.squared
    ## Collect data to anotab returned as the 'anova' object
    anotab <-  list()
    ## Step forward and continue as long as R2.adj improves and R2.adj
    ## remains below R2.adj < R2.all
    R2.previous <- R2.0
    repeat {
        if (trace) {
            cat("Step: R2.adj=", R2.previous, "\n")
            cat(pasteCall(formula(object)), "\n")
        }
        adds <- add.scope(object, scope)
        ## Nothing to add and we're done: break
        if (length(adds) == 0)
            break
        R2.adds <- numeric(length(adds))
        adds <- paste("+", adds)
        names(R2.adds) <- adds
        ## Loop over add scope
        for (trm in seq_along(R2.adds)) {
            fla <- paste(". ~ .", names(R2.adds[trm]))
            R2.tmp <- suppressMessages(RsquareAdj(update(object, fla),
                                 permutations = R2permutations, ...)$adj.r.squared)
            if (!length(R2.tmp) || is.na(R2.tmp))
                R2.tmp <- 0
            R2.adds[trm] <- R2.tmp
        }
        best <- which.max(R2.adds)
        if (trace) {
            out <- sort(c("<All variables>" = R2.all, "<model>" = R2.previous,
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
            (!R2scope || R2scope && R2.adds[best] <= R2.all)) {
            ## Second criterion: added variable is significant
            tst <- suppressMessages(
                add1(object, scope = adds[best], test="permu",
                     permutations = permutations,
                     alpha = Pin, trace = FALSE, ...))
            if (trace) {
                print(tst[-1,])
                cat("\n")
            }
            if (tst[,"Pr(>F)"][2] <= Pin) {
                fla <- paste("~  .", names(R2.adds[best]))
                object <-  update(object, fla)
            } else
                break
        } else {
            break
        }
        R2.previous <- RsquareAdj(object,
                                  permutations = R2permutations,
                                  ...)$adj.r.squared
            anotab <- rbind(anotab,
                            cbind("R2.adj" = R2.previous, tst[2,]))
    }
    if (NROW(anotab)) {
        if (R2scope)
            anotab <- rbind(anotab, "<All variables>" = c(R2.all, rep(NA, 4)))
        class(anotab) <- c("anova", class(anotab))
        object$anova <- anotab
    }
    object
}
