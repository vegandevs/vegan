`density.anosim` <-
    function(x, ...)
{
    out <- density(x$perm, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}
