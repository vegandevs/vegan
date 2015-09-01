### Casts a vegan spantree object into single linkage dendrogram of
### class hclust. The non-trivial items in "hclust" object are a
### 'merge' matrix for fusions of points and/or clusters, a 'height'
### vector which gives the heights of each fusion, and an 'order'
### vector that gives the order of leaves in the plotted
### dendrogram. The 'height's are only sorted spantree segment
### distances, but for 'merge' we need to establish cluster
### memberships, and for 'order' we must traverse the tree.

`as.hclust.spantree` <-
    function(x, ...)
{
    ## Order by the lengths of spanning tree links
    o <- order(x$dist)
    npoints <- x$n
    if(npoints < 2)
        stop("needs at least two points")
    ## Ordered indices of dads and kids
    dad <- (2:npoints)[o]
    kid <- x$kid[o]
    ## merge matrix of hclust has negative index when a single point
    ## is added to a tree and a positive index when a group is joined
    ## to a tree, and the group is numbered by the level it was
    ## formed.
    labs <- -seq_len(npoints)
    merge <- matrix(0, nrow=npoints-1, ncol=2)
    for(i in 1:nrow(merge)) {
        merge[i, ] <- c(labs[dad[i]], labs[kid[i]])
        ## update labs for the current group and its kids
        labs[labs %in% labs[c(dad[i], kid[i])]] <- i
    }

    order <- hclustMergeOrder(merge)
      
    out <- list(merge = merge, height = x$dist[o], order = order,
                labels = x$labels, method = "spantree", call =
                match.call())
    class(out) <- "hclust"
    out
}

### Internal vegan function to get the 'order' from a merge matrix of
### an hclust tree

`hclustMergeOrder` <-
    function(merge)
{
    ## Get order of leaves with recursive search from the root
    order <- numeric(nrow(merge)+1)
    ind <- 0
    ## "<<-" updates data only within hclustMergeOrder, but outside
    ## the visit() function.
    visit <- function(i, j) {
        if (merge[i,j] < 0) {
            ind <<- ind+1
            order[ind] <<- -merge[i,j]
        } else {
            visit(merge[i,j], 1)
            visit(merge[i,j], 2)
        }
    }
    visit(nrow(merge), 1)
    visit(nrow(merge), 2)
    return(order)
}

### Reorder an hclust tree. Basic R provides reorder.dendrogram, but
### this functoin works with 'hclust' objects, and also differs in
### implementation. We use either weighted mean, min or max or
### sum. The dendrogram is always ordered in ascending order, so that
### with max the left kid always has lower value. So with 'max' the
### largest value is smaller in leftmost group. The choice 'sum'
### hardly makes sense, but it is the default in
### reorder.dendrogram. The ordering with 'mean' differs from
### reorder.dendrogram which uses unweighted means, but here we weight
### means by group sizes so that the mean of an internal node is the
### mean of its leaves.

`reorder.hclust` <-
    function(x, wts,
             agglo.FUN = c("mean", "min", "max", "sum", "uwmean"),
             ...)
{
    agglo.FUN <- match.arg(agglo.FUN)
    merge <- x$merge
    nlev <- nrow(merge)
    stats <- numeric(nlev)
    counts <- numeric(nlev)
    pair <- numeric(2)
    pairw <- numeric(2)
    ## Go through merge, order each level and update the statistic.
    for(i in 1:nlev) {
        for(j in 1:2) {
            if (merge[i,j] < 0) {
                pair[j] <- wts[-merge[i,j]]
                pairw[j] <- 1
            } else {
                pair[j] <- stats[merge[i,j]]
                pairw[j] <- counts[merge[i,j]]
            }
        }
        ## reorder
        merge[i,] <- merge[i, order(pair)]
        ## statistic for this merge level
        stats[i] <-
            switch(agglo.FUN,
                   "mean" = weighted.mean(pair, pairw),
                   "min" = min(pair),
                   "max" = max(pair),
                   "sum" = sum(pair),
                   "uwmean" = mean(pair))
        counts[i] <- sum(pairw)
    }
    ## Get the 'order' of the reordered dendrogram
    order <- hclustMergeOrder(merge)
    x$merge <- merge
    x$order <- order
    x$value <- stats
    x
}

### Trivial function to reverse the order of an hclust tree (why this
### is not in base R?)

`rev.hclust` <-
    function(x)
{
    x$order <- rev(x$order)
    x
}

### Get coordinates for internal or terminal nodes (leaves) that would
### be used in plot.hclust

`scores.hclust` <-
    function(x, display = "internal", ...)
{
    extnam <- c("leaves", "terminal")
    intnam <- c("internal")
    display <- match.arg(display, c(extnam, intnam))
    ## Terminal nodes (leaves): plot.hclust scales x-axis for n points
    ## as 1..n. The y-value is the 'height' where the terminal node
    ## was fused to the tree.
    if(display %in% extnam) {
        merge <- x$merge
        y <- numeric(nrow(merge) + 1)
        for(i in 1:nrow(merge))
            for(j in 1:2)
                if(merge[i,j] < 0)
                    y[-merge[i,j]] <- x$height[i]
        xx <- order(x$order)
        xy <- cbind(`x` = xx, `height` = y)
    } else {
        ## Internal nodes are given in the order they were fused which
        ## also is the order of 'height'
        xx <- reorder(x, order(x$order), agglo.FUN = "uwmean")$value
        xy <- cbind(`x`= xx, `height` = x$height)
    }
    xy
}
