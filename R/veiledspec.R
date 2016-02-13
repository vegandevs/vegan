"veiledspec" <-
    function(x, ...)
{
    if (!inherits(x, "prestonfit"))
        x <- prestonfit(x)
    S.obs <- sum(x$freq)
    p <- x$coefficients
    S.tot <- p["S0"]*p["width"]*sqrt(2*pi)
    out <- c(S.tot, S.obs, S.tot - S.obs)
    names(out) <- c("Extrapolated","Observed","Veiled")
    out
}
