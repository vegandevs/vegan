`density.permutest.cca` <-
    function(x, ...)
{
    out <- density(x$F.perm, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}
