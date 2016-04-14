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

## cca result: no RsquareAdj
RsquareAdj.cca <-
    function(x, nperm, print_progress=TRUE, ...)
{
    r2 <- x$CCA$tot.chi/x$tot.chi
    n <- nrow(x$CCA$u)
    rand_r2 <- rep(NA, nperm)
    Y_string <- as.character(x$terms[[2]])
    Y <- eval(parse(text=Y_string))
    x$call[2] <- sub(Y_string, 'Y_rand', x$call[2])
    for (i in 1:nperm) {
        Y_rand <- Y[sample(n), ]
        cca_rand <- eval(parse(text=paste(x$call[1], '(',x$call[2], 
                                         ', data=', x$call[3], ')', 
                                         sep='')))
        rand_r2[i] <- cca_rand$CCA$tot.chi / x$tot.chi
        if (print_progress) {
            if (i %% 100 == 0)  
                print(paste('perm:', i, 'adj.r.squared:', 
                            round(1 - ((1 - r2) / (1 - mean(rand_r2, na.rm=T))),3)))
        }
    }
    # Eq 5 Peres-Neto et al. 2006
    radj <- 1 - ((1 - r2) / (1 - mean(rand_r2)))
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
