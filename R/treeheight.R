`treeheight` <-
    function(tree)
{
    if (inherits(tree, "spantree"))
        return(sum(tree$dist))
    tree <- as.hclust(tree)
    ## nodes should start from 0 -- if there are negative heights,
    ## tree is too pathological to be measured.
    if (any(tree$height < 0))
        stop("negative heights: tree cannot be measured")
    ## can be done really fast if there are no reversals, but we need
    ## to traverse the tree with reversals
    if (is.unsorted(tree$height)) { # slow
        h <- tree$height
        m <- tree$merge
        height <- 0
        for (i in 1:nrow(m)) {
            for (j in 1:2) {
                if (m[i,j] < 0)
                    height <- height + h[i]
                else
                    height <- height + abs(h[i] - h[m[i,j]])
            }
        }
        height
    }
    else    # fast
        sum(tree$height) + max(tree$height)
}

