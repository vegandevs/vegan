`permutest.betadisper` <- function(x, pairwise = FALSE,
                                   control = permControl(nperm = 999), ...)
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
    if(!inherits(x, "betadisper"))
        stop("Only for class \"betadisper\"")
    ## will issue error if only a single group
    mod.aov <- anova(x)
    nobs <- length(x$distances)
    mod <- lm(x$distances ~ x$group)
    mod.Q <- mod$qr
    p <- mod.Q$rank
    resids <- qr.resid(mod.Q, x$distances)
    res <- numeric(length = control$nperm + 1)
    res[1] <- summary(mod)$fstatistic[1]
    ## pairwise comparisons
    if(pairwise) {
        ## unique pairings
        combin <- combn(levels(x$group), 2)
        n.pairs <- ncol(combin)
        t.stats <- matrix(0, ncol = n.pairs, nrow = control$nperm + 1)
        t.stats[1,] <- apply(combn(levels(x$group), 2), 2, function(z) {
            t.statistic(x$distances[x$group == z[1]],
                        x$distances[x$group == z[2]])})
    }
    for(i in seq(along = res[-1])) {
        perm <- shuffle(nobs, control = control)
        perm.resid <- resids[perm]
        f <- qr.fitted(mod.Q, perm.resid)
        mss <- sum((f - mean(f))^2)
        r <- qr.resid(mod.Q, perm.resid)
        rss <- sum(r^2)
        rdf <- nobs - p
        resvar <- rss / rdf
        res[i+1] <- (mss / (p - 1)) / resvar
        ## pairwise comparisons
        if(pairwise) {
            for(j in seq_len(n.pairs)) {
                grp1 <- x$distance[perm][x$group == combin[1, j]]
                grp2 <- x$distance[perm][x$group == combin[2, j]]
                t.stats[i+1, j] <- t.statistic(grp1, grp2)
            }
        }
    }
    pval <- sum(res >= res[1]) / length(res)
    if(pairwise) {
        df <- apply(combin, 2, function(z) {
            length(x$distances[x$group == z[1]]) +
                length(x$distance[x$group == z[2]]) - 2})
        pairwise <- list(observed = 2 * pt(-abs(t.stats[1,]), df),
                         permuted = apply(t.stats, 2,
                         function(z) sum(abs(z) >= abs(z[1]))/length(z)))
        names(pairwise$observed) <- names(pairwise$permuted) <-
            apply(combin, 2, paste, collapse = "-")
    } else {
        pairwise <- NULL
    }
    retval <- cbind(mod.aov[, 1:4], c(control$nperm, NA), c(pval, NA))
    dimnames(retval) <- list(c("Groups", "Residuals"),
                             c("Df", "Sum Sq", "Mean Sq", "F", "N.Perm",
                               "Pr(>F)"))
    retval <- list(tab = retval, pairwise = pairwise,
                   groups = levels(x$group), control = control)
    class(retval) <- "permutest.betadisper"
    retval
}
