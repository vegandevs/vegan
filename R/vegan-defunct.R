## as.mlm: deprecated in 2.5-0, defunct in 2.6-0

`as.mlm` <-
    function(x)
{
    .Defunct("see ?hatvalues.cca for new alternatives")
    if (is.null(x$CCA))
        stop("'as.mlm' can be used only for constrained ordination")
    UseMethod("as.mlm")
}
`as.mlm.cca` <-
    function (x)
{
    w <- x$rowsum
    WA <- x$CCA$wa
    X <- qr.X(x$CCA$QR)
    ## shall use weighted regression: deweight X
    X <- (1/sqrt(w)) * X
    X <- as.data.frame(X)
    lm(WA ~ ., data = X, weights = w)
}

`as.mlm.rda` <-
    function (x)
{
    X <- as.data.frame(qr.X(x$CCA$QR))
    WA <- x$CCA$wa
    lm(WA ~ . , data = X)
}
### commsimulator was deprecated in 2.4-0, defunct in 2.6-0

"commsimulator" <-
function (x, method, thin = 1)
{
    .Defunct("simulate(nullmodel(x, method))", package="vegan")
    method <- match.arg(method,
                        c("r0","r1","r2","r00","c0","swap", "tswap",
                          "backtrack", "quasiswap"))
    if (method == "r0")
        method <- "r0_old"
    x <- as.matrix(x)
    out <- simulate(nullmodel(x, method), nsim = 1, thin = thin)
    out <- out[,,1]
    attributes(out) <- attributes(x)
    out
}
