# Get Species or Site Scores from an Ordination

Function to access either species or site scores for specified axes in
some ordination methods. The `scores` function is generic in vegan, and
vegan ordination functions have their own `scores` functions that are
documented separately with the method (see e.g.
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`scores.metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md),
[`scores.decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)).
This help file documents the default `scores` method that is only used
for non-vegan ordination objects.

## Usage

``` r
# Default S3 method
scores(x, choices,
    display=c("sites", "species", "both"), tidy = FALSE, ...)
```

## Arguments

- x:

  An ordination result.

- choices:

  Ordination axes. If missing, default method returns all axes.

- display:

  Partial match to access scores for `"sites"` or `"species"` of for
  `"both"`.

- tidy:

  Return `"both"` scores in data frame that is compatible with
  [ggplot2](https://CRAN.R-project.org/package=ggplot2), with variable
  `score` labelling the scores as `"sites"` or `"species"`.

- ...:

  Other parameters (unused).

## Details

Function `scores` is a generic method in vegan. Several vegan functions
have their own `scores` methods with their own defaults and with some
new arguments. This help page describes only the default method. For
other methods, see, e.g.,
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`scores.rda`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`scores.decorana`](https://vegandevs.github.io/vegan/reference/decorana.md).

All vegan ordination functions should have a `scores` method which
should be used to extract the scores instead of directly accessing them.
Scaling and transformation of scores should also happen in the `scores`
function. If the `scores` function is available, the results can be
plotted using
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md),
[`ordixyplot`](https://vegandevs.github.io/vegan/reference/ordixyplot.md)
etc., and the ordination results can be compared in
[`procrustes`](https://vegandevs.github.io/vegan/reference/procrustes.md)
analysis.

The `scores.default` function is used to extract scores from non-vegan
ordination results. Many standard ordination methods of libraries do not
have a specific `class`, and no specific method can be written for them.
However, `scores.default` guesses where some commonly used functions
keep their site scores and possible species scores.

If `x` is a matrix, `scores.default` returns the chosen columns of that
matrix, ignoring whether species or sites were requested (do not regard
this as a bug but as a feature, please). Currently the function seems to
work at least for [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html),
[`prcomp`](https://rdrr.io/r/stats/prcomp.html),
[`princomp`](https://rdrr.io/r/stats/princomp.html) and some ade4
objects. It may work in other cases or fail mysteriously.

## Value

The function returns a matrix of scores if one type is requested, or a
named list of matrices if `display = "both"`, or a
[ggplot2](https://CRAN.R-project.org/package=ggplot2) compatible data
frame if `tidy = TRUE`.

## Author

Jari Oksanen

## See also

Specific `scores` functions include (but are not limited to)
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`scores.rda`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`scores.decorana`](https://vegandevs.github.io/vegan/reference/decorana.md),
[`scores.envfit`](https://vegandevs.github.io/vegan/reference/envfit.md),
[`scores.metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md),
[`scores.monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md)
and
[`scores.pcnm`](https://vegandevs.github.io/vegan/reference/pcnm.md).
These have somewhat different interface –
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
in particular – but all work with keywords `display="sites"` and return
a matrix. However, they may also return a list of matrices, and some
other `scores` methods will have quite different arguments.

## Examples

``` r
data(varespec)
vare.pca <- prcomp(varespec)
scores(vare.pca, choices=c(1,2))
#>            PC1         PC2
#> 18 -10.7847878  18.7094315
#> 15 -27.8036826 -11.7414745
#> 24 -25.6919559 -14.5399684
#> 27 -31.7820166 -31.2216800
#> 23 -19.6315869  -2.5541193
#> 19  -0.2413294 -11.4974077
#> 22 -26.6771373 -12.3140897
#> 16 -21.9230366   0.4449159
#> 28 -39.6083051 -41.8877392
#> 13  -4.0664328  20.4191153
#> 14 -18.4416245   5.4406988
#> 20 -17.3999191   2.3653380
#> 25 -25.1673547 -13.2508067
#> 7  -11.4065430  41.7356300
#> 5   -8.4243752  45.3805255
#> 6   -2.0759474  36.9311222
#> 3   39.8617580   8.0590041
#> 4   13.1065901  12.8377217
#> 2   57.6827011  -4.8983565
#> 9   63.3138332 -22.4481549
#> 12  44.1073111 -10.1653935
#> 10  64.9418975 -16.7633564
#> 11  11.5313633   3.9720890
#> 21  -3.4194194  -3.0130455
```
