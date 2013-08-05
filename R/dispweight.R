`dispweight` <-
    function(comm, group, nperm = 1000)
{
    if (missing(group))
        group <- rep(1, nrow(comm))
    # number of replicates per group
    nrep <- tabulate(group) 
    # workhorse
    dfun <- function(comm, group, nperm, nrep) {
        ## Calc Dispersion
        # group means
        means <-  tapply(comm, group, mean)
        # omit groups with mean == 0
        if(any(means == 0)) {
            comm <- comm[means[group] != 0]
            group <- group[means[group] != 0, drop = TRUE]
            nrep <- nrep[means != 0]
            means <- means[means != 0]
        }
        # group variances
        vars <-  tapply(comm, group, var)
        # dispersion
        d <- vars / means
        # average dispersion 
        d_hat <- sum(d * (nrep - 1)) / sum(nrep - 1)
        
        ## Test
        # original chisq values
        chi_o <- sum((comm - means[group])^2 / means[group])
        
        # permutations
        # calculate chisq
        pfun <- function(comm, group){
            means <-  tapply(comm, group, mean)
            if(any(means == 0)) {
                comm <- comm[means[group] != 0]
                group <- group[means[group] != 0, drop = TRUE]
                means <- means[means != 0]
            }
            chi <- sum((comm - means[group])^2 / means[group])
            return(chi)
        }
        # sum of individuals per group
        sums <- tapply(comm, group, sum)
        # realocate randomly individuals to replications
        perms <- vector('list', length(sums))
        for(i in seq_along(sums)) {
            perms[[i]] <- rmultinom(n = nperm, size = sums[i], prob = rep(1, nrep[i]) / nrep[i])
        }
        perms <- t(do.call(rbind, perms))
        chi_p <- apply(perms, 1, pfun, sort(group))
        p <- (sum(chi_p >= chi_o) + 1) / (nperm + 1)
        out <- list(D = d_hat, p = p, weights = ifelse(p < 0.05, 1/d_hat, 1))
        return(out)
    }
    # apply workhorse to every species
    out <- apply(comm, 2, dfun, group, nperm, nrep)
    
    # format output
    weights <-  unlist(sapply(out, '[', 3))
    out = sweep(comm, MARGIN = 2, weights, `*`)
    attr(out, "D") <- unlist(sapply(out, '[', 1)) 
    attr(out, "p") <- unlist(sapply(out, '[', 2))
    attr(out, "weights") <-  weights
    attr(out, "permutations") <- nperm
    class(out) <- c("dispweight", class(out))
    return(out)
}
