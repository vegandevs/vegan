`dispweight` <-
    function(comm, groups, nsimul = 1000)
{
    ## only applicable for counts
    if (!identical(all.equal(comm, round(comm)), TRUE))
        stop("function needs counts (integers)")
    ## no groups
    if (missing(groups))
        groups <- rep(1, nrow(comm))
    # number of replicates per group
    nrep <- tabulate(groups) 
    # workhorse
    dfun <- function(comm, groups, nsimul, nrep) {
        ## Calc Dispersion
        # group means
        means <-  tapply(comm, groups, mean)
        # omit groups with mean == 0
        if(any(means == 0)) {
            comm <- comm[means[groups] != 0]
            groups <- groups[means[groups] != 0, drop = TRUE]
            nrep <- nrep[means != 0]
            means <- means[means != 0]
        }
        # group variances
        vars <-  tapply(comm, groups, var)
        # dispersion
        d <- vars / means
        # average dispersion 
        d_hat <- sum(d * (nrep - 1)) / sum(nrep - 1)
        
        ## Test
        # original chisq values
        chi_o <- sum((comm - means[groups])^2 / means[groups])
        
        # permutations
        # calculate chisq
        pfun <- function(comm, groups){
            means <-  tapply(comm, groups, mean)
            if(any(means == 0)) {
                comm <- comm[means[groups] != 0]
                groups <- groups[means[groups] != 0, drop = TRUE]
                means <- means[means != 0]
            }
            chi <- sum((comm - means[groups])^2 / means[groups])
            return(chi)
        }
        # sum of individuals per group
        sums <- tapply(comm, groups, sum)
        # realocate randomly individuals to replications
        perms <- vector('list', length(sums))
        for(i in seq_along(sums)) {
            perms[[i]] <- rmultinom(n = nsimul, size = sums[i], prob = rep(1, nrep[i]) / nrep[i])
        }
        perms <- t(do.call(rbind, perms))
        chi_p <- apply(perms, 1, pfun, sort(groups))
        p <- (sum(chi_p >= chi_o) + 1) / (nsimul + 1)
        out <- list(D = d_hat, p = p, weights = ifelse(p < 0.05, 1/d_hat, 1))
        return(out)
    }
    # apply workhorse to every species
    out <- apply(comm, 2, dfun, groups, nsimul, nrep)
    
    # format output
    weights <-  unlist(sapply(out, '[', 3))
    out = sweep(comm, MARGIN = 2, weights, `*`)
    attr(out, "D") <- unlist(sapply(out, '[', 1)) 
    attr(out, "p") <- unlist(sapply(out, '[', 2))
    attr(out, "weights") <-  weights
    attr(out, "nsimul") <- nsimul
    class(out) <- c("dispweight", class(out))
    return(out)
}
