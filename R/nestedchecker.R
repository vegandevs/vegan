"nestedchecker" <-
function(comm)
{
    cb <- sum(designdist(comm, "(A-J)*(B-J)", "binary"))
    sppairs <- ncol(comm)*(ncol(comm)-1)/2
    out <- list("C.score" = cb/sppairs, statistic = cb)
    class(out) <- "nestedchecker"
    out
}

