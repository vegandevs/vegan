### The depths of nodes in a 'spantree' object: The nodes are either
### leaves with one link, or internal nodes with >1 links. The leaves
### are removed recursively from the tree and at each step the depth
### is increased with one.
`spandepth` <-
    function (x) 
{
    if (!inherits(x, "spantree"))
        stop("'x' must be 'spantree' result")
    kid <- c(NA, x$kid)
    par <- p <- seq_along(kid)
    par[1] <- NA
    ## Isolated nodes in disconnected tree have depth 0, other nodes
    ## start from depth 1
    intree <- p %in% kid | !is.na(kid) 
    depth <- numeric(length(par))
    depth[intree] <- 1
    if (!is.null(x$labels))
        names(depth) <- x$labels
    while(any(intree)) {
        ## Node is internal (intree) if it is both a parent and a kid
        ## and kid is in the tree or it is kid to two or more parents
        intree <- (p %in% intersect(kid[intree], par[intree]) &
                   p %in% p[intree][kid[intree] %in% p[intree]] |
                   p %in% kid[intree][duplicated(kid[intree])]) 
        depth[intree] <- depth[intree] + 1
    }
    depth
}

