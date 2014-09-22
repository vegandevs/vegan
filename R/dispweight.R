`dispweight` <-
    function(comm, groups, nsimul = 999, nullmodel = "c0_ind",
             plimit = 0.05)
{
    ## no groups?
    if (missing(groups))
        groups <- rep(1, nrow(comm))
    ## Remove empty levels of 'groups' or this fails cryptically (and
    ## take care 'groups' is a factor)
    groups <- factor(groups)
    ## Statistic is the sum of squared differences by 'groups'
    means <- apply(comm, 2, function(x) tapply(x, groups, mean))
    ## handle 1-level factors: all sites belong to the same 'groups'
    if (is.null(dim(means)))
        means <- matrix(means, nrow=1, ncol = length(means),
                        dimnames = list(levels(groups), names(means)))
    ## expand to matrix of species means
    fitted <- means[groups,]
    dhat <- colSums((comm - fitted)^2/fitted, na.rm = TRUE)
    ## Get df for non-zero blocks of species. Completely ignoring
    ## all-zero blocks for species sounds strange, but was done in the
    ## original paper, and we follow here. However, this was not done
    ## for significance tests, and only concerns 'D' and 'weights'.
    nreps <- table(groups)
    div <- colSums(sweep(means > 0, 1, nreps - 1, "*"))
    ## "significance" of overdispersion is assessed from Chi-square
    ## evaluated separately for each species. This means fixing only
    ## marginal totals for species but letting row marginals vary
    ## freely, unlike in standard Chi-square where both margins are
    ## fixed. In vegan this is achieved by nullmodel 'c0_ind'. Instead
    ## of one overall simulation, nullmodel is generated separately
    ## for each of 'groups'
    chisq <- function(x) {
        fitted <- colMeans(x)
        colSums(sweep(x, 2, fitted)^2, na.rm = TRUE) / fitted
    }
    simulated <- matrix(0, nrow = ncol(comm), ncol = nsimul)
    for (lev in levels(groups)) {
        nm <- nullmodel(comm[groups == lev,], nullmodel)
        if (nm$commsim$binary)
            stop("'binary' nullmodel cannot be used")
        tmp <- apply(simulate(nm, nsimul), 3, chisq)
        ok <- !is.na(tmp)
        simulated[ok] <- simulated[ok] + tmp[ok] 
    }
    ## p value based on raw dhat, then we divide
    p <- (rowSums(dhat <= simulated) + 1) / (nsimul + 1)
    dhat <- dhat/div
    weights <- ifelse(p <= plimit, 1/dhat, 1)
    comm <- sweep(comm, 2, weights, "*")
    attr(comm, "D") <- dhat
    attr(comm, "df") <- div
    attr(comm, "p") <- p
    attr(comm, "weights") <-  weights
    attr(comm, "nsimul") <- nsimul
    attr(comm, "nullmodel") <- nullmodel
    class(comm) <- c("dispweight", class(comm))
    comm
}
