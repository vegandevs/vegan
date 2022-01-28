### unbiased Simpson index, Hurlbert (1971) "nonconcept" paper, eq. 5,
### but implemented here with rarefy (because I'm lazy and just re-use
### work already done).

`simpson.unb` <-
    function(x, inverse = FALSE)
{
    d <- rarefy(x, 2) - 1
    ## alternatively use directly the Hurlbert equation
    ## n <- rowSums(x)
    ## d <- rowSums(x/n*(n-x)/(n-1))
    if (inverse)
        d <- 1/(1-d)
    attr(d, "Subsample") <- NULL
    d
}
