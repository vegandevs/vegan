`nesteddisc` <-
    function(comm, niter = 200)
{
    ## The original discrepancy method orders columns by frequencies,
    ## but does not consider ties. The current function tries to order
    ## tied values to minimize the discrepancy either by complete
    ## enumeration or with a larger number of ties using simulated
    ## annealing.  NALL: max no. of tied items for NALL! complete
    ## enumeration

    ## starting values and CONSTANTS
    NALL <- 7
    allperm <- factorial(NALL)
    ties <- FALSE
    trace <- FALSE
    ## Code
    comm <- ifelse(comm > 0, 1, 0)
    cs <- colSums(comm)
    k <- rev(order(cs))
    ## initial order
    cs <- cs[k]
    comm <- comm[, k]
    ## run lengths: numbers of tied values
    le <- rle(cs)$lengths
    cle <- c(0, cumsum(le))
    x <- seq(along=cs)
    ## Range of row sums: only swaps between these have an effect
    rs <- range(rowSums(comm))
    ## Function to evaluate discrepancy
    FUN <- function(x) sum(comm[col(comm)[,x] <= rowSums(comm)] == 0)
    Ad <- FUN(x)
    ## Go through all le-items and permute ties. Function shuffleSet
    ## is in permute package, and its minperm argument triggers
    ## complete enumeration for no. of ties <= NALL!.
    for (i in seq_along(le)) {
        if (le[i] > 1) {
            take <- x
            idx <- seq_len(le[i]) + cle[i]
            ## Can swaps influence discrepancy?
            if (idx[1] > rs[2] || idx[le[i]] < rs[1])
                next
            perm <- shuffleSet(le[i], niter, control = how(minperm = allperm),
                               quietly = TRUE)
            ## maxperm is a double -- needs EPS -0.5
            if ((attr(perm, "control")$maxperm - 0.5) > niter)
                ties <- TRUE
            perm <- matrix(perm, ncol = le[i]) + cle[i]
            vals <- sapply(1:nrow(perm), function(j) {
                take[idx] <- perm[j,]
                FUN(take)
            })
            jmin <- which.min(vals)
            if (vals[jmin] < Ad) {
                x[idx] <- perm[jmin,]
                Ad <- vals[jmin]
            }
        }
    }
    out <- list(statistic=Ad, ties = ties, order = k[x])
    names(out$statistic) <- "discrepancy"
    class(out) <- "nesteddisc"
    out
}

