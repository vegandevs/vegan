## as.mlm: deprecated in 2.5-0, defunct in 2.6-0

`as.mlm` <-
    function(x)
{
    .Defunct("see ?hatvalues.cca for new alternatives")
}

`as.mlm.cca` <-
    function(x)
{
    .Defunct("see ?hatvalues.cca for new alternatives")
}

### commsimulator was deprecated in 2.4-0, defunct in 2.6-0

"commsimulator" <-
function (x, method, thin = 1)
{
    .Defunct("simulate(nullmodel(x, method))", package="vegan")
}

### deprecated in 2.2-0, but forgotten and never exported from the
### NAMESPACE. Make finally defunct for 2.6-0.

"permuted.index" <-
    function (n, strata)
{
    .Defunct("use permute package (shuffle or shuffleSet)")
}

### deprecated in 2.5-3, declared defunct and removed in 2.6-2 (but
### actually was not removed)
`humpfit` <-
    function (mass, spno, family = poisson, start)
{
    .Defunct("use natto::humpfit() from https://github.com/jarioksa/natto/")
}

### deprecated in 2.6-2 with alternative toCoda, make defunct in 2.6-5
`as.mcmc.oecosimu` <-
    function(x)
{
    ## Deprecated in favour of toCoda: using as an S3 method would
    ## need importFrom(coda, as.mcmc) and that would add dependence on
    ## coda
    .Defunct("toCoda", package = "vegan")
}

`as.mcmc.permat` <- as.mcmc.oecosimu
