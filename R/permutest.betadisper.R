`permutest.betadisper` <- function(x, pairwise = FALSE,
                                   permutations = 999,
                                   parallel = getOption("mc.cores"), ...)
{
    EPS <- sqrt(.Machine$double.eps) # for P-value comparisons
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
            dim(idx) <- c(1, length(idx))
        }
        R <- nrow(idx)
        Fperm <- matrix(nrow = R, ncol = 1)
        if (pairwise) {                 # set up object to hold t stats
            Tperm <- matrix(ncol = n.pairs, nrow = R)
            Jseq <- seq_len(n.pairs)
        }
        rdf <- nobs - p                 # residual degrees of freedom
        ## iterate
        for (i in seq_len(R)) {         # iterate
            take <- idx[i, ]                # current permutation from set
            p.resid <- resids[take]         # permute residuals
            f <- qr.fitted(mod.Q, p.resid)  # create new data
            mss <- sum((f - mean(f))^2)
            r <- qr.resid(mod.Q, p.resid)
            rss <- sum(r^2)
            resvar <- rss / rdf
            Fperm[i, ] <- (mss / (p - 1)) / resvar

            ## pairwise tests
            if(pairwise) {
                for(j in Jseq) {
                    grp1 <- x$distance[take][group == combin[1, j]]
                    grp2 <- x$distance[take][group == combin[2, j]]
                    Tperm[i, j] <- t.statistic(grp1, grp2)
                }
            }
        }

        ## bind on pairwise stats if any
        if (pairwise) {
            Fperm <- cbind(Fperm, Tperm)
        }
        Fperm
    }

    if(!inherits(x, "betadisper"))
        stop("only for class \"betadisper\"")

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

    ## pairwise comparisons
    if(pairwise) {
        combin <- combn(levels(x$group), 2) # unique pairings
        n.pairs <- ncol(combin)
    }

    ## Parallel processing of permutations
    if (is.null(parallel)) {
        parallel <- 1
    }
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1L) {
        if (.Platform$OS.type == "unix" && !hasClus) {
            Pstats <- do.call("rbind",
                           mclapply(seq_len(nperm),
                                    function(x) permFun(permutations[x, , drop = FALSE]),
                                    mc.cores = parallel))
        } else {
            ## if hasClus, don't set up and top a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
            }

            Pstats <- parApply(parallel, permutations, 1,
                               function(x) permFun(x))
            if (is.null(dim(Pstats))) {
                Pstats <- matrix(Pstats) # one-column matrix
            } else {
                Pstats <- t(Pstats) # transpose statistics to columns
            }

            if (!hasClus) {
                stopCluster(parallel)
            }
        }
    } else {
        Pstats <- permFun(permutations)
    }

    ## Process results
    F0 <- summary(mod)$fstatistic[1]
    Fstats <- Pstats[, 1, drop = TRUE]    # allow empty dim to be dropped
    statistic <- F0
    names(statistic) <- "Overall (F)"

    ## pairwise comparisons
    if(pairwise) {
        T0 <- apply(combn(levels(group), 2), 2, function(z) {
            t.statistic(x$distances[group == z[1]],
                        x$distances[group == z[2]])})
        Tstats <- Pstats[, -1, drop = FALSE]
        statistic <- c(statistic, T0)
    }

    ## compute permutation p-value
    pval <- (sum(Fstats >= F0 - EPS) + 1) / (length(Fstats) + 1)

    if(pairwise) {
        df <- apply(combin, 2, function(z) {
            length(x$distances[group == z[1]]) +
                length(x$distance[group == z[2]]) - 2})
        pairp <- (colSums(sweep(abs(Tstats), 2, abs(T0), '>=')) + 1) /
            (NROW(Tstats) + 1)
        pairp <- list(observed = 2 * pt(-abs(T0), df),
                         permuted = pairp)
        tnames <- apply(combin, 2, paste, collapse = "-")
        names(pairp$observed) <- names(pairp$permuted) <- tnames
        names(statistic)[-1] <- paste(tnames, "(t)")
    } else {
        pairp <- NULL
    }

    retval <- cbind(mod.aov[, 1:4], c(nperm, NA), c(pval, NA))
    dimnames(retval) <- list(c("Groups", "Residuals"),
                             c("Df", "Sum Sq", "Mean Sq", "F", "N.Perm",
                               "Pr(>F)"))
    retval <- list(tab = retval,
                   pairwise = pairp,
                   groups = levels(group),
                   statistic = statistic,
                   perm = if (pairwise) {
                       structure(cbind(Fstats, Tstats), dimnames = list(NULL, names(statistic)))
                   } else {
                       structure(Fstats, names = names(statistic))
                   },
                   control = attr(permutations, "control"))
    class(retval) <- "permutest.betadisper"
    retval
}
