### rewritten commsimulator

"commsimulator" <-
function (x, method, thin = 1) 
{
    ## Do not yet warn on deprecation to allow smooth transition
    .Deprecated("simulate(nullmodel(x, method))", package="vegan")
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
