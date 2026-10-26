`nestedchecker` <-
    function(comm)
{
    cb <- sum(designdist(comm, "(A-J)*(B-J)", "binary"))
    sppairs <- ncol(comm)*(ncol(comm)-1)/2
    out <- list("C.score" = cb/sppairs, statistic = as.integer(cb))
    names(out$statistic) <- "checkerboards"
    class(out) <- "nestedchecker"
    out
}

