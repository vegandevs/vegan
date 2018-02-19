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

## cca result: RsquareAdj is calculated similarly as in
## varpart(). This is similar "semipartial" model as for rda() and
## found as a difference of R2-adj values of combined model with
## constraints + conditions and only conditions.
`RsquareAdj.cca` <-
    function (x, permutations = 1000, ...)
{
    r2 <- x$CCA$tot.chi / x$tot.chi
    if (is.null(x$pCCA)) {
        p <- permutest(x, permutations, ...)
        radj <- 1 - ((1 - r2) / (1 - mean(p$num / x$tot.chi)))
    } else {
        p <- getPermuteMatrix(permutations, nobs(x))
        Y <- ordiYbar(x, "initial")
        r2tot <- (x$pCCA$tot.chi + x$CCA$tot.chi) / x$tot.chi
        r2null <- mean(sapply(seq_len(nrow(p)), function(i)
            sum(qr.fitted(x$CCA$QR, Y[p[i,],])^2)))
        r2tot <- 1 - ((1-r2tot)/(1-r2null/x$tot.chi))
        r2p <- x$pCCA$tot.chi / x$tot.chi
        r2null <- mean(sapply(seq_len(nrow(p)), function(i)
            sum(qr.fitted(x$pCCA$QR, Y[p[i,],])^2)))
        r2p <- 1 - ((1-r2p)/(1-r2null/x$tot.chi))
        radj <- r2tot - r2p
    }
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
