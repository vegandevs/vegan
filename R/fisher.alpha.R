`fisher.alpha` <-
    function (x, MARGIN = 1, ...)
{
    x <- as.matrix(x, rownames.force = TRUE)
    if(ncol(x) == 1)
        x <- t(x)
    sol <- apply(x, MARGIN, fisherfit)
    unlist(lapply(sol, function(x) x$estimate))
}
