`treeheight` <-
    function(tree)
{
    if (inherits(tree, "spantree"))
        return(sum(tree$dist))
    tree <- as.hclust(tree)
    sum(tree$height) + max(tree$height)
}

