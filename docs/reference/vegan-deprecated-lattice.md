# Deprecated lattice Functions in vegan

Lattice functions are going to be deprecated and removed from vegan.
They are replaced with ggplot2 functions in CRAN package
[ggvegan](https://CRAN.R-project.org/package=ggvegan).

## Usage

``` r
ordicloud(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi3d", prepanel = "prepanel.ordi3d", ...)
ordisplom(x, data = NULL, formula = NULL, display = "sites", choices = 1:3, 
    panel = "panel.ordi", type = "p", ...)
ordiresids(x, kind = c("residuals", "scale", "qqmath"),
    residuals = "working", type = c("p", "smooth", "g"),
    formula, ...)
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
```

## Arguments

- x:

  Input object.

- kind:

  The type of plot: residuals or absolute values of residuals against
  fitted values, or quantile plot of residuals with
  [`qqmath`](https://rdrr.io/pkg/lattice/man/qqmath.html).

- residuals:

  The type of residuals with choices `"working"`, `"response"`,
  `"standardized"` and `"studentized"`.

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

- panel, prepanel:

  The name of the panel or prepanel function.

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

Trellis (or lattice) functions were added to vegan mostly in 2008 to
2009. In that time they were the only alternative of the kind, but now
there are better, more versatile and more user-friendly alternatives,
mainly in ggplot2. CRAN package ggvegan provides modern alternatives to
most lattice functions in vegan. The lattice functions in vegan will be
deprecated as soon ggvegan provides a ggplot2 alternative. The
deprecated functions will be defunct in the next major release of vegan.

The following functions are currently deprecated:

- `ordicloud` was transferred to
  [vegan3d](https://CRAN.R-project.org/package=vegan3d) as
  `ordilattice3d`.

- `ordisplom` design is bad and deficient. If you want to have something
  similar, write your own code.

- `ordixyplot`: use `autoplot` or `ordiggplot` functions in ggvegan.

- `plot` functions for
  [`poolaccum`](https://vegandevs.github.io/vegan/reference/specpool.md),
  [`renyi`](https://vegandevs.github.io/vegan/reference/renyi.md) and
  [`renyiaccum`](https://vegandevs.github.io/vegan/reference/renyi.md):
  use `autoplot` in ggvegan.

- `permulattice`: use `autoplot` in ggvegan.

- `ordiresids` is not very useful, but you can directly access
  ordination results with
  [`fitted.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
  [`residuals.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
  [`rstandard.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
  [`rstudent.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md)
  and other functions that were not available in `ordiresids`.
