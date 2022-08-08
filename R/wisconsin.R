`wisconsin` <-
    function(x)
{
    x <- decostand(x, "max", 2)
    mx <- attr(x, "parameters")$max
    x <- decostand(x, "total", 1)
    attr(x, "parameters")$max <- mx
    attr(x, "decostand") <- "wisconsin"
    x
}
