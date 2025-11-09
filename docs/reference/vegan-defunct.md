# Defunct Functions in Package vegan

The functions or variables listed here are no longer part of vegan as
they are no longer needed.

## Usage

``` r
## defunct in vegan 2.8-0
ordicloud(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi3d", prepanel = "prepanel.ordi3d", ...)
ordisplom(x, data = NULL, formula = NULL, display = "sites", choices = 1:3, 
    panel = "panel.ordi", type = "p", ...)
ordiresids(x, kind = c("residuals", "scale", "qqmath"),
    residuals = "working", type = c("p", "smooth", "g"),
    formula, ...)

## defunct in vegan 2.7-0
adonis(formula, data, permutations = 999, method = "bray",
    strata = NULL, contr.unordered = "contr.sum",
    contr.ordered = "contr.poly", parallel = getOption("mc.cores"), ...)
orditkplot(...)

## defunct in vegan 2.6-0
as.mlm(x)
humpfit(mass, spno, family = poisson, start)
vegandocs(doc = c("NEWS", "ONEWS", "FAQ-vegan", "intro-vegan",
    "diversity-vegan", "decision-vegan", "partitioning", "permutations"))

## defunct in vegan 2.5-0
commsimulator(x, method, thin=1)

## defunct in vegan 2.4-0
# S3 method for class 'adonis'
density(x, ...)
# S3 method for class 'vegandensity'
plot(x, main = NULL, xlab = NULL, ylab = "Density", 
   type = "l", zero.line = TRUE, obs.line = TRUE, ...)
# S3 method for class 'adonis'
densityplot(x, data, xlab = "Null", ...)

## defunct in vegan 2.2-0
metaMDSrotate(object, vec, na.rm = FALSE, ...)

## defunct in vegan 2.0-0
getNumObs(object, ...)
permuted.index2(n, control = permControl())
```

## Details

Lattice function`ordicloud` was moved to
[vegan3d](https://CRAN.R-project.org/package=vegan3d) as
`ordilattice3d`. To substitute `ordiresids`, use
[`influence.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md)
to extract data that you can use to design your own graphics.
`ordisplom` was unsatisfactory and there is no replacement

[`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.md)
replaces `adonis` with extended functionality and completely new
internal design. The shared arguments of `adonis` are similar as in
[`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.md), but
arguments `contr.unordered` and `contr.ordered` can set the contrasts
within `adonis`.

`orditkplot` was moved to CRAN package
[vegan3d](https://CRAN.R-project.org/package=vegan3d) version 1.3-0.
Install vegan3d and use the function in the old way.

`as.mlm` function is replaced with a set functions that can find the
same statistics directly from the ordination result object: see
[`influence.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md).

Function `humpfit` was transferred to the natto package and is still
available from <https://github.com/jarioksa/natto/>.

R functions [`news`](https://rdrr.io/r/utils/news.html) should be used
to read vegan NEWS (`news(package = "vegan")`), and
[`browseVignettes`](https://rdrr.io/r/utils/browseVignettes.html) is a
better tool for reading vignettes than `vegandocs`.

Function `commsimulator` is replaced with
[`make.commsim`](https://vegandevs.github.io/vegan/reference/commsim.md)
which defines the Null models, and functions
[`nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
and
[`simulate.nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
that check the input data and generate the Null model communities.

The deprecated `density` and `densityplot` methods are replaced with
similar methods for
[`permustats`](https://vegandevs.github.io/vegan/reference/permustats.md).
The
[`permustats`](https://vegandevs.github.io/vegan/reference/permustats.md)
offers more powerful analysis tools for permutations, including
[`summary.permustats`](https://vegandevs.github.io/vegan/reference/permustats.md)
giving \\z\\ values (a.k.a. standardized effect sizes, SES), and Q-Q
plots
([`qqnorm.permustats`](https://vegandevs.github.io/vegan/reference/permustats.md),
[`qqmath.permustats`](https://vegandevs.github.io/vegan/reference/permustats.md)).

Function `metaMDSrotate` is replaced with
[`MDSrotate`](https://vegandevs.github.io/vegan/reference/MDSrotate.md)
which can handle
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md)
results in addition to
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md).

The permutation functions were moved to the permute package, and they
are documented there. The permute package replaces `permuted.index` and
`permuted.index2` with
[`shuffle`](https://rdrr.io/pkg/permute/man/shuffle.html) and
`getNumObs` with its specific
[`nobs-methods`](https://rdrr.io/pkg/permute/man/nobs.html).

## See also

[`Defunct`](https://rdrr.io/r/base/Defunct.html),
[`vegan-deprecated`](https://vegandevs.github.io/vegan/reference/vegan-deprecated.md)
