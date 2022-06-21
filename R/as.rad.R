`as.rad` <-
    function(x)
{
    if (inherits(x, "rad"))
        return(x)
    ## recursive call for several observations
    if (isTRUE(nrow(x) > 1)) {
        comm <- apply(x, 1, as.rad)
        class(comm) <- "rad.frame"
        return(comm)
    }
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

## do not print 'index' attribute

`print.rad` <-
    function(x, ...)
{
    print(as.table(x), ...)
    invisible(x)
}
