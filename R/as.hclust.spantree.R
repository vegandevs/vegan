### Casts a vegan spantree object into single linkage dendrogram of
### class hclust. The non-trivial items in "hclust" object are a
### 'merge' matrix for fusions of points and/or clusters, a 'height'
### vector which gives the heights of each fusion, and an 'order'
### vector that gives the order of leaves in the plotted
### dendrogram. The 'height's are only sorted spantree segment
### distances, but for 'merge' we need to establish cluster
### memberships, and for 'order' we must traverse the tree. The
### plot.hclust() function seems to require that the left kid is
### always more compact (a single point or fused earlier than the
### right kid).

`as.hclust.spantree` <-
    function(x, ...)
{
    ## Order by the lengths of spanning tree links
    o <- order(x$dist)
    npoints <- length(o) + 1
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
        ## add items of labs, keep tighter cluster on the left
        items <- c(labs[dad[i]], labs[kid[i]])
        if (items[1] > 0 || items[2] > 0)
            items <- sort(items)
        else
            items <- rev(sort(items))
        merge[i, ] <- items
        ## update labs for the current group and its kids
        labs[labs %in% labs[c(dad[i], kid[i])]] <- i
    }
    ## Get order of leaves with recursive search from the root
    visited <- matrix(FALSE, nrow = nrow(merge), ncol=ncol(merge))
    order <- numeric(npoints)
    ind <- 0
    ## "<<-" updates data only within this function, but outside the
    ## visit() function.
    visit <- function(i, j) {
        if(visited[i,j])
            return(NULL)
        else {
            visited[i,j] <<- TRUE
        }
        if (merge[i,j] < 0) {
            ind <<- ind+1
            order[ind] <<- -merge[i,j]
            if (j == 1)
                visit(i, 2)
        } else {
            visit(merge[i,j], 1)
            visit(merge[i,j], 2)
        }
    }
    visit(nrow(merge), 1)
    visit(nrow(merge), 2)
    
    out <- list(merge = merge, height = x$dist[o], order = order,
                labels = x$labels, method = "spantree", call =
                match.call())
    class(out) <- "hclust"
    out
}
