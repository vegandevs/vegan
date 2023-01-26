`bstick.decorana` <-
    function(n, ...)
{
    tot.chi <- n$totchi
    ## assume full rank input
    n.comp <- min(nrow(n$rproj), nrow(n$cproj)) - 1
    res <- bstick.default(n.comp, tot.chi, ...)
    ## only four axes in decorana
    res <- res[1:4]
    names(res) <- names(n$evals)
    res
}

