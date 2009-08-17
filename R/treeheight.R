`treeheight` <-
    function(tree)
{
    if (inherits(tree, "spantree"))
        return(sum(tree$dist))
    tree <- as.hclust(tree)
    m <- tree$merge
    h <- tree$height
    height <- 0
    for (i in 1:nrow(m)) {
        for (j in 1:2) {
            if (m[i,j] < 0)
                height <- height + h[i]
            else
                height <- height + h[i] - h[m[i,j]] 
        }
    }
    height
}

