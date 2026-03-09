# Deprecated Functions in vegan package

These functions are provided for compatibility with older versions of
vegan only, and may be defunct as soon as the next release.

## Usage

``` r
## Following lattice functions are deprecated. Use ggplot2 functions
## in the CRAN package ggvegan instead

ordixyplot(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi", aspect = "iso", envfit,
    type = c("p", "biplot"), ...)
# S3 method for class 'poolaccum'
plot(x, alpha = 0.05, type = c("l","g"), ...)
# S3 method for class 'renyi'
plot(x, ...)
# S3 method for class 'renyiaccum'
plot(x, what = c("Collector", "mean", "Qnt 0.025",
    "Qnt 0.975"),
     type = "l", ...)
permulattice(x, plot = c("densityplot", "qqmath"), observed = TRUE,
    axislab = "Permutations", ...)
# S3 method for class 'permustats'
densityplot(x, data, observed = TRUE,
    xlab = "Permutations", ...)
# S3 method for class 'permustats'
qqmath(x, data, observed = TRUE, sd.scale = FALSE,
    ylab = "Permutations", ...)

## use toCoda instead
as.mcmc.oecosimu(x)
as.mcmc.permat(x)
```

## Arguments

- x:

  input object

- data:

  Optional data to amend ordination results. The ordination results are
  found from `x`, but you may give here data for other variables needed
  in plots. Typically these are environmental data.

- formula:

  Formula to define the plots. A default formula will be used if this is
  omitted. The ordination axes must be called by the same names as in
  the ordination results (and these names vary among methods).

- display:

  The kind of scores: an argument passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md).

- choices:

  The axes selected: an argument passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md).

- panel:

  The name of the panel function.

- aspect:

  The aspect of the plot (passed to the lattice function).

- envfit:

  Result of
  [`envfit`](https://vegandevs.github.io/vegan/reference/envfit.md)
  function displayed in `ordixyplot`. Please note that this needs same
  `choices` as `ordixyplot`.

- type:

  The type of plot. This knows the same alternatives as
  [`panel.xyplot`](https://rdrr.io/pkg/lattice/man/panel.xyplot.html).
  In addition `ordixyplot` has alternatives `"biplot"`, `"arrows"` and
  `"polygon"`. The first displays fitted vectors and factor centroids of
  `envfit`, or in constrained ordination, the biplot arrows and factor
  centroids if `envfit` is not given. The second (`type = "arrows"`) is
  a trellis variant of
  [`ordiarrows`](https://vegandevs.github.io/vegan/reference/ordiarrows.md)
  and draws arrows by `groups`. The line parameters are controlled by
  [`trellis.par.set`](https://rdrr.io/pkg/lattice/man/trellis.par.get.html)
  for `superpose.line`, and the user can set `length`, `angle` and
  `ends` parameters of
  [`panel.arrows`](https://rdrr.io/pkg/lattice/man/llines.html). The
  last one (`type = "polygon"`) draws a polygon enclosing all points in
  a panel over a polygon enclosing all points in the data. The overall
  polygon is controlled by Trellis parameters
  [`trellis.par.set`](https://rdrr.io/pkg/lattice/man/trellis.par.get.html)
  `plot.polygon` and `superpose.polygon`.

- what:

  Items to be plotted.

- plot:

  Use lattice function
  [`densityplot`](https://rdrr.io/pkg/lattice/man/histogram.html) or
  [`qqmath`](https://rdrr.io/pkg/lattice/man/qqmath.html).

- xlab, ylab, axislab:

  Label for the axis displaying permutation values.

- observed:

  Add observed statistic among permutations.

- sd.scale:

  Scale permutations to unit standard deviation and observed statistic
  to standardized effect size.

- alpha:

  Level of quantiles shown. This proportion will be left outside
  symmetric limits.

- ...:

  Arguments passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md)
  methods or lattice functions.

## Details

`as.mcmc` functions were replaced with
[`toCoda`](https://vegandevs.github.io/vegan/reference/oecosimu.md).

Better alternatives to lattice functions are provided by ggplot2
functions in the [ggvegan](https://CRAN.R-project.org/package=ggvegan)
package. `ggvegan::autoplot` functions are most similar to these
deprecated functions.

## See also

[`Deprecated`](https://rdrr.io/r/base/Deprecated.html)
