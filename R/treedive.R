`treedive` <-
    function(comm, tree, match.force = FALSE)
{
    if (!inherits(tree, c("hclust", "spantree")))
        stop("'clus' must be an 'hclust' or 'spantree' result object")
    m <- as.matrix(cophenetic(tree))
    ## Check tree/comm match by names
    if (match.force || ncol(comm) != ncol(m)) {
        if (match.force)
            message("Forced matching of 'tree' labels and 'comm' names")
        else
            message("Dimensions do not match between 'comm' and 'tree'")
        fnd <- colnames(comm) %in% tree$labels
        if (!all(fnd)) {
            warning("not all names of 'comm' found in 'tree'")
            comm <- comm[, fnd]
        }
        fnd <- tree$labels %in% colnames(comm)
        if (!all(fnd))
            warning("not all names of 'tree' found in 'comm'")
        comm <- comm[, tree$labels[fnd]]
        if (length(unique(tree$labels)) != length(tree$labels))
            stop("names not unique in 'tree': match wrong")
        if (length(unique(colnames(comm))) != ncol(comm))
            stop("names not unique in 'comm': match wrong")
    }
    ## Repeat for sites
    div <- numeric(nrow(comm))
    for (i in 1:nrow(comm)) {
        k <- comm[i,] > 0
        nit <- sum(k)
        ## Trivial cases of zero or one species
        if (nit==0)
            div[i] <- NA
        else if (nit==1)
            div[i] <- 0
        else {
            d <- as.dist(m[k,k])
            cl <- update(tree, d = d)
            div[i] <- treeheight(cl)
        }
    }
    names(div) <- rownames(comm)
    div
}
