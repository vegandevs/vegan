`bstick.cca` <-
    function(n, ...)
{
    if(!inherits(n, c("rda", "cca")))
        stop("'n' not of class \"cca\" or \"rda\"")
    if(!is.null(n$CCA) && n$CCA$rank > 0)
        stop("'bstick' only for unconstrained models")
    ## No idea how to define bstick for dbrda or capscale with
    ## negative eigenvalues
    if (inherits(n, c("dbrda", "capscale")) &&
        (!is.null(n$CA$imaginary.u) || !is.null(n$CA$imaginary.u.eig)))
        stop(gettextf("'bstick' cannot be used for '%s' with negative eigenvalues",
                      class(n)[1]))
    ## need to select appropriate total inertia
    tot.chi <- n$CA$tot.chi
    n.comp <- n$CA$rank
    res <- bstick.default(n.comp, tot.chi, ...)
    names(res) <- names(n$CA$eig)
    res
}
