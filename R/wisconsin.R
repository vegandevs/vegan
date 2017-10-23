`wisconsin` <-
    function(x)
{
    x <- decostand(x, "max", 2)
    decostand(x, "tot", 1)
}
