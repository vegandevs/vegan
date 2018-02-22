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
        p <- x$pCCA$QR$rank
        radj <- RsquareAdj(R2 + R2p, n, m + p) - RsquareAdj(R2p, n, p)
    }
    list(r.squared = R2, adj.r.squared = radj)
}

## cca result: R2_adj calculation for partial models differs from
## varpart(), where it is assessed as difference of R2_adj values for
## conditions and constraints. Here we use difference of unadjusted R2
## values, and adjust that difference.
`RsquareAdj.cca` <-
    function (x, permutations = 1000, ...)
{
    r2 <- x$CCA$tot.chi / x$tot.chi
    ## model='direct' makes a difference in pCCA: it guarantees that
    ## dependent data are permuted in the same way in the partial and
    ## constrained components so that we can take their difference for
    ## the correction term.
    p <- permutest(x, permutations, model = "direct")
    radj <- 1 - ((1 - r2) / (1 - mean(p$num) / x$tot.chi))
    list(r.squared = r2, adj.r.squared = radj)
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
