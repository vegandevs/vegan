`specnumber` <-
    function(x, groups, MARGIN = 1)
{
    if (!missing(groups)) {
        if (length(groups) == 1)
            groups <- rep(groups, nrow(x))
        x <- aggregate(x, list(groups), max)
        rownames(x) <- x[,1]
        x <- x[,-1]
    }
    if (length(dim(x)) > 1)
        apply(x > 0, MARGIN, sum)
    else
        sum(x > 0)
}
