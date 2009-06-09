`RsquareAdj` <-
    function(x, ...)
{
    UseMethod("RsquareAdj")
}


`RsquareAdj.default` <-
    function(x, n, m, ...)
{
    if (m >= (n-1))
        NA
    else
        1 - (1-x)*(n-1)/(n-m-1)
}

## Use this with rda() results
`RsquareAdj.rda` <-
    function(x, ...)
{
    R2 <- x$CCA$tot.chi/x$tot.chi
    m <- x$CCA$rank
    n <- nrow(x$CCA$u)
    if (is.null(x$pCCA))
        radj <- RsquareAdj(R2, n, m)
    else
        radj <- NA
    list(r.squared = R2, adj.r.squared = radj)
}

## cca result: no RsquareAdj
RsquareAdj.cca <-
    function(x, ...)
{
    R2 <- x$CCA$tot.chi/x$tot.chi
    m <- x$CCA$rank
    n <- nrow(x$CCA$u)
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
