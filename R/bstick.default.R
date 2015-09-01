`bstick.default` <-
    function(n, tot.var = 1, ...)
{
    res <- rev(cumsum(tot.var/n:1)/n)
    names(res) <- paste("Stick", seq(len=n), sep="")
    res
}
