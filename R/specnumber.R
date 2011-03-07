`specnumber` <-
    function(x, MARGIN = 1)
{
    if (length(dim(x)) > 1)
        apply(x > 0, MARGIN, sum)
    else
        sum(x > 0)
}
