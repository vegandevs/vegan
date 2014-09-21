`permutest.betadisper` <- function(x, pairwise = FALSE,
                                   permutations = how(nperm = 999),
                                   parallel = getOption("mc.cores"), ...)
{
    t.statistic <- function(x, y) {
        m <- length(x)
        n <- length(y)
        xbar <- mean(x)
        ybar <- mean(y)
        xvar <- var(x)
        yvar <- var(y)
        pooled <- sqrt(((m-1)*xvar + (n-1)*yvar) / (m+n-2))
        (xbar - ybar) / (pooled * sqrt(1/m + 1/n))
    }

    permFun <- function(idx) {
        if (!is.matrix(idx)) {
            dim(i) <- c(1, length(idx))
        }
        R <- nrow(idx)
        Fperm <- numeric(length = R)
        ## iterate
        for (i in seq_len(R)) {
            take <- idx[i, ]                # current permutation from set
            p.resid <- resids[take]         # permute residuals
            f <- qr.fitted(mod.Q, p.resid) # create new data
            mss <- sum((f - mean(f))^2)
            r <- qr.resid(mod.Q, p.resid)
            rss <- sum(r^2)
            rdf <- nobs - p
            resvar <- rss / rdf
            Fperm[i] <- (mss / (p - 1)) / resvar
        }
        Fperm
    }

    if(!inherits(x, "betadisper"))
        stop("Only for class \"betadisper\"")

    ## will issue error if only a single group
    mod.aov <- anova(x)
    nobs <- length(x$distances) ## number of observations
    mod <- lm(x$distances ~ x$group)
    mod.Q <- mod$qr
    p <- mod.Q$rank
    resids <- qr.resid(mod.Q, x$distances)

    ## extract groups
    group <- x$group
    ## permutations is either a single number, a how() structure or a
    ## permutation matrix
    permutations <- getPermuteMatrix(permutations, nobs)
    nperm <- nrow(permutations)

    ## Parallel processing of permutations
    if (is.null(parallel)) {
        parallel <- 1
    }
    hasClus <- inherits(parallel, "cluster")
    if ((hasClus || parallel > 1L) && requireNamespace("parallel")) {
        if (.Platform$OS.type == "unix" && !hasClus) {
            Fstats <- do.call("c",
                           mclapply(seq_len(nperm),
                                    function(x) permFun(permutations[x, , drop = FALSE]),
                                    mc.cores = parallel))
        } else {
            ## if hasClus, don't set up and top a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
            }
            Fstats <- parRapply(parallel, permutations, function(x) permFun(x))
            if (!hasClus) {
                stopCluster(parallel)
            }
        }
    } else {
        Fstats <- permFun(permutations)
    }
    F0 <- summary(mod)$fstatistic[1]
    Fstats <- round(Fstats, 12)
    F0 <- round(F0, 12)

    ## set-up objects to hold permuted results
    ## res <- numeric(length = nperm + 1)
    ## res[1] <- summary(mod)$fstatistic[1]
    ## pairwise comparisons
    if(pairwise) {
        ## unique pairings
        combin <- combn(levels(x$group), 2)
        n.pairs <- ncol(combin)
        t.stats <- matrix(0, ncol = n.pairs, nrow = nperm + 1)
        t.stats[1,] <- apply(combn(levels(group), 2), 2, function(z) {
            t.statistic(x$distances[group == z[1]],
                        x$distances[group == z[2]])})
    }

    ## begin loop over shuffleSet perms
    for(i in seq_len(nperm)) {
        perm <- permutations[i,] ## take current permutation from set
        ## perm.resid <- resids[perm] ## permute residuals
        ## f <- qr.fitted(mod.Q, perm.resid) ## create new data
        ## mss <- sum((f - mean(f))^2)
        ## r <- qr.resid(mod.Q, perm.resid)
        ## rss <- sum(r^2)
        ## rdf <- nobs - p
        ## resvar <- rss / rdf
        ## res[i+1] <- (mss / (p - 1)) / resvar

        ## pairwise comparisons
        if(pairwise) {
            for(j in seq_len(n.pairs)) {
                grp1 <- x$distance[perm][group == combin[1, j]]
                grp2 <- x$distance[perm][group == combin[2, j]]
                t.stats[i+1, j] <- t.statistic(grp1, grp2)
            }
        }
    }

    ## compute permutation p-value
    ## pval <- sum(res >= res[1]) / length(res)
    pval <- (sum(Fstats >= F0) + 1) / (length(Fstats) + 1)

    if(pairwise) {
        df <- apply(combin, 2, function(z) {
            length(x$distances[group == z[1]]) +
                length(x$distance[group == z[2]]) - 2})
        pairwise <- list(observed = 2 * pt(-abs(t.stats[1,]), df),
                         permuted = apply(t.stats, 2,
                         function(z) sum(abs(z) >= abs(z[1]))/length(z)))
        names(pairwise$observed) <- names(pairwise$permuted) <-
            apply(combin, 2, paste, collapse = "-")
    } else {
        pairwise <- NULL
    }

    retval <- cbind(mod.aov[, 1:4], c(nperm + 1, NA), c(pval, NA))
    dimnames(retval) <- list(c("Groups", "Residuals"),
                             c("Df", "Sum Sq", "Mean Sq", "F", "N.Perm",
                               "Pr(>F)"))
    retval <- list(tab = retval, pairwise = pairwise,
                   groups = levels(group),
                   control = attr(permutations, "control"))
    class(retval) <- "permutest.betadisper"
    retval
}
