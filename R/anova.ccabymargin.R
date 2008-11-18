`anova.ccabymargin` <-
    function(object, step=100, scope, ...)
{
    if (!missing(scope) && is.character(scope))
        trms <- scope
    else
        trms <- drop.scope(object, scope)
    alltrms <- labels(terms(object$terminfo))
    keep <- trms %in% alltrms
    trms <- trms[keep]
    ntrms <- length(trms)
    bigperm <- 0
    for (i in 1:ntrms) {
        fla <- formula(object)
        ## Put all trms except current into Condition() and update
        ## formula
        if (length(alltrms) > 1) {
            keeptrms <- alltrms[!(alltrms==trms[i])]
            updfla <- paste("Condition(",paste(keeptrms, collapse="+"), ")")
            fla <- update(fla, paste(". ~ . + ", updfla))
        }
        tmp <- update(object, fla)
        tmp <- anova(tmp, step=step, ...)
        ## Start every permutation from the same seed, but get the
        ## seed of the longest simulation and reset the RNG to that
        ## state when exiting the function
        if (tmp[1,"N.Perm"] > bigperm) {
            bigperm <- tmp[1, "N.Perm"]
            bigseed <- get(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)
        }
        if (i == 1) {
            seed <- attr(tmp, "Random.seed")
            sol <- tmp
        }
        else {
            sol <- rbind(sol[1:(i-1),], as.matrix(tmp[1,]), sol[i,])
            assign(".Random.seed", seed, envir = .GlobalEnv)
        }
    }
    ## Put RNG at the end of the longest simulation
    assign(".Random.seed", bigseed, envir = .GlobalEnv)
    rownames(sol)[1:ntrms] <- trms
    head <- attr(sol, "heading")
    head[1] <- paste(head[1], "Marginal effects of terms\n", sep="")
    head[2] <- paste("Model:", c(object$call))
    attr(sol, "heading") <- head
    sol
}
