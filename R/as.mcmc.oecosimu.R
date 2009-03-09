`as.mcmc.oecosimu` <-
    function(x) 
{
    x <- as.ts(x)
    mcpar <- attr(x, "tsp")
    mcpar[3] <- round(1/mcpar[3])
    attr(x, "mcpar") <- mcpar
    class(x) <- c("mcmc", class(x))
    x
}
