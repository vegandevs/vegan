"specnumber" <-
    function(x, MARGIN = 1)
{
    apply(x > 0, MARGIN, sum)
}
