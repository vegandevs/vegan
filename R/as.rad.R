"as.rad" <-
    function(x)
{
    if (inherits(x, "rad"))
        return(x)
    take <- x > 0
    nm <- names(x)
    comm <- x[take]
    names(comm) <- nm[take]
    comm <- rev(sort(comm))
    class(comm) <- "rad"
    comm
}
