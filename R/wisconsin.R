`wisconsin` <-
    function(x, na.rm = FALSE)
{
    x <- decostand(x, "max", 2, na.rm = na.rm)
    mx <- attr(x, "parameters")$max
    x <- decostand(x, "total", 1, na.rm = na.rm)
    attr(x, "parameters")$max <- mx
    attr(x, "decostand") <- "wisconsin"
    x
}
