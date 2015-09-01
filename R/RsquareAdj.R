`RsquareAdj` <-
    function(x, ...)
{
    UseMethod("RsquareAdj")
}


`RsquareAdj.default` <-
    function(x, n, m, ...)
{
    r2 <- 1 - (1-x)*(n-1)/(n-m-1)
    if (any(na <- m >= n-1))
        r2[na] <- NA
    r2
}

## Use this with rda() results
`RsquareAdj.rda` <-
    function(x, ...)
{
    R2 <- x$CCA$tot.chi/x$tot.chi
    m <- x$CCA$qrank
    n <- nrow(x$CCA$u)
    if (is.null(x$pCCA)) {
        radj <- RsquareAdj(R2, n, m)
    } else {
        ## Partial model: same adjusted R2 as for component [a] in two
        ## source varpart model
        R2p <- x$pCCA$tot.chi/x$tot.chi
        p <- x$pCCA$rank
        radj <- RsquareAdj(R2 + R2p, n, m + p) - RsquareAdj(R2p, n, p)
    }
    list(r.squared = R2, adj.r.squared = radj)
}

## dbRDA: Euclidean style distances with no imaginary component can be
## handled as rda, but I have no idea how to handle objects with
## imaginary inertia.

`RsquareAdj.capscale` <-
    function(x, ...)
{
    if (!is.null(x$CA$imaginary.chi))
        list(r.squared = NA, adj.r.squared = NA)
    else
        NextMethod("RsquareAdj", x, ...)
}

## cca result: no RsquareAdj
RsquareAdj.cca <-
    function(x, ...)
{
    R2 <- x$CCA$tot.chi/x$tot.chi
    radj <- NA
    list(r.squared = R2, adj.r.squared = radj)
}

## Linear model: take the result from the summary
RsquareAdj.lm <-
  function(x, ...)
{
    summary(x)[c("r.squared", "adj.r.squared")]
}

## Generalized linear model: R2-adj only with Gaussian model
RsquareAdj.glm <-
    function(x, ...)
{
    if (family(x)$family == "gaussian")
        summary.lm(x)[c("r.squared", "adj.r.squared")]
    else
        NA
}
