## deprecated in 2.7-0, defunct in 2.8-0
`ordicloud` <-
    function(x, data = NULL, formula, display = "sites", choices=1:3,
             panel = "panel.ordi3d",
             prepanel = "prepanel.ordi3d", ...)
{
    .Defunct("vegan3d::ordilattice3d", package="vegan")
}

`ordiresids` <-
    function(x, kind = c("residuals", "scale", "qqmath"),
             residuals = "working",
             type = c("p", "smooth", "g"), formula, ...)
{
    .Defunct("influence.cca to extract residuals etc.", package="vegan")
}

`ordisplom` <-
    function(x, data=NULL, formula = NULL,  display = "sites", choices = 1:3,
             panel = "panel.ordi", type = "p", ...)
{
    .Defunct(package="vegan")
}

## adonis: replaced by adonis2
## adonis2 added in 2.4-0
## adonis message on future deprecation in 2.6-2, deprecated in 2.6-6
##        defunct in 2.7-1
`adonis` <-
    function(formula, data=NULL, permutations=999, method="bray", strata=NULL,
             contr.unordered="contr.sum", contr.ordered="contr.poly",
             parallel = getOption("mc.cores"), ...)
{
    .Defunct("adonis2", package="vegan")
}

## as.mlm: deprecated in 2.5-0, defunct in 2.6-0

`as.mlm` <-
    function(x)
{
    .Defunct("influence.cca", package = "vegan")
}

`as.mlm.cca` <-
    function(x)
{
    .Defunct("influence.cca", package = "vegan")
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
    .Defunct("permute package (shuffle or shuffleSet)")
}

### deprecated in 2.5-3, declared defunct and removed in 2.6-2 (but
### actually was not removed)
`humpfit` <-
    function (mass, spno, family = poisson, start)
{
    .Defunct("natto::humpfit() from https://github.com/jarioksa/natto/")
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
