### Multiple-site dissimilarity indices (Sorensen & Jaccard) and their
### decomposition into "turnover" and "nestedness" following Baselga
### (Global Ecology & Biogeography 19, 134-143; 2010). Implemented as
### nestedness functions and directly usable in oecosimu().

`nestedbetasor` <-
    function(comm)
{
    beta <- betadiver(comm, method = NA)
    b <- beta$b
    c <- beta$c
    diffbc <- sum(abs(b-c))
    sumbc <- sum(b+c)
    bmin <- sum(pmin(b, c))
    a <- sum(comm) - sum(colSums(comm) > 0)
    simpson <- bmin/(bmin + a)
    nest <- a/(bmin + a) * diffbc/(2*a + sumbc)
    sorensen <- sumbc/(2*a + sumbc)
    c(turnover = simpson, nestedness = nest, sorensen = sorensen)
}

`nestedbetajac` <-
    function(comm)
{
    beta <- betadiver(comm, method = NA)
    b <- beta$b
    c <- beta$c
    diffbc <- sum(abs(b-c))
    sumbc <- sum(b+c)
    bmin <- sum(pmin(b, c))
    a <- sum(comm) - sum(colSums(comm) > 0)
    simpson <- 2*bmin/(2*bmin + a)
    nest <- a/(2*bmin + a) * diffbc/(a + sumbc)
    jaccard <- sumbc/(a + sumbc)
    c(turnover = simpson, nestedness = nest, jaccard = jaccard)
}
