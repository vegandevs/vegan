`bstick.princomp` <-
    function(n, ...)
{
    if(!inherits(n, "princomp"))
        stop("'n' not of class \"princomp\"")
    tot.chi <- sum(n$sdev^2)
    n.comp <- length(n$sdev)
    res <- bstick.default(n.comp, tot.chi, ...)
    names(res) <- dimnames(n$loadings)[[2]]
    res
}
