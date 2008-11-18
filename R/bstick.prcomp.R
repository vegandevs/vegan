`bstick.prcomp` <-
    function(n, ...)
{
    if(!inherits(n, "prcomp"))
        stop("'n' not of class \"prcomp\"")
    tot.chi <- sum(n$sdev^2)
    n.comp <- length(n$sdev)
    res <- bstick.default(n.comp, tot.chi, ...)
    names(res) <- dimnames(n$rotation)[[2]]
    res
}
