`as.rad` <-
    function(x)
{
    if (inherits(x, "rad"))
        return(x)
    take <- x > 0
    nm <- names(x)
    comm <- x[take]
    names(comm) <- nm[take]
    comm <- sort(comm, decreasing = TRUE, index.return = TRUE)
    ## ordered index of included taxa
    index <- which(take)[comm$ix]
    comm <- comm$x
    attr(comm, "index") <- index
    class(comm) <- "rad"
    comm
}
