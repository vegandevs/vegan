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

## rda. Function varpart() estimates adjusted R2 for partial model
## rda(Y, X1, X2) as a difference of adjusted R2 values of models
## rda(Y, cbind(X1,X2)) and rda(Y, X2), and we did the same here from
## vegan release 2.0-3 to 2.4. However, I am not quite convinced that
## this is correct, and therefore I disable adjusted R2 for partial
## models while I study the issue (Jari Oksanen, 23/2/2018).  If you
## change handling of partial models, remember to update ordiR2step!
`RsquareAdj.rda` <-
    function(x, ...)
{
    R2 <- x$CCA$tot.chi/x$tot.chi
    m <- x$CCA$qrank
    n <- nrow(x$CCA$u)
    if (is.null(x$pCCA)) {
        radj <- RsquareAdj(R2, n, m)
    } else {
        radj <- NA
    }
    list(r.squared = R2, adj.r.squared = radj)
}

## cca: for pCCA, see comment for rda. Earlier we used the same
## equation both for pCCA and CCA which was inconsistent with rda. git
## commit 9c25e6fdb4 has a version that was consistent with rda in
## vegan 2.0-3 to 2.4-6.
`RsquareAdj.cca` <-
    function (x, permutations = 1000, ...)
{
    r2 <- x$CCA$tot.chi / x$tot.chi
    if (is.null(x$pCCA) || x$pCCA$rank == 0) {
        p <- permutest(x, permutations, ...)
        radj <- 1 - ((1 - r2) / (1 - mean(p$num) / x$tot.chi))
    } else {
        radj <- NA
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
