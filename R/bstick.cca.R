`bstick.cca` <-
    function(n, ...)
{
    if(!inherits(n, c("rda", "cca")))
        stop("'n' not of class \"cca\" or \"rda\"")
    if(!is.null(n$CCA))
        stop("'bstick' only for unconstrained models.")
    ## need to select appropriate total invertia
    tot.chi <- n$CA$tot.chi
    n.comp <- n$CA$rank
    res <- bstick.default(n.comp, tot.chi, ...)
    names(res) <- names(n$CA$eig)
    res
}
