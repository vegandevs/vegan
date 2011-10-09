`nesteddisc` <-
    function(comm)
{
    ## The original discrepancy method orders columns by frequencies,
    ## but does not consider ties. The current function tries to order
    ## tied values to minimize the discrepancy either by complete
    ## enumeration or with a larger number of ties using simulated
    ## annealing.  NALL: max no. of tied items for NALL! complete
    ## enumeration

    ## starting values and CONSTANTS
    NALL <- 7
    NITER <- 200
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
    ## Go through all le-items and permute ties
    for (i in 1:length(le)) {
        if (le[i] > 1) {
            take <- x
            idx <- (1:le[i]) + cle[i]
            ## Can swaps influence discrepancy?
            if (idx[1] > rs[2] || idx[le[i]] < rs[1])
                next
            Ad <- FUN(x)
            ## Complete enumeration if no. of tied value <= NALL
            if (le[i] <= NALL) {
                perm <- matrix(allPerms(le[i]), ncol=le[i]) + cle[i]
                ## Take at maximum NITER cases from complete enumeration
                if (nrow(perm) >= NITER) {
                    perm <- perm[sample.int(nrow(perm), NITER),]
                    ties <- TRUE
                }
            }
            ## No complete enumeration, but a sample and potentially
            ## duplicated orders
            else {
                ties <- TRUE
                perm <- matrix(0, nrow=NITER, ncol=le[i])
                for (j in 1:NITER)
                    perm[j,] <- permuted.index(le[i])
                perm <- perm + cle[i]
            }
            for (j in 1:nrow(perm)) {
                take[idx] <- perm[j,]
                val <- FUN(take)
                if (val < Ad) {
                    x <- take
                    Ad <- val
                    if (trace)
                        cat(Ad, ":", perm[j,], "\n")
                }
            }
        }
    }
    out <- list(statistic=Ad, ties = ties, order = k[x])
    class(out) <- "nesteddisc"
    out
}

