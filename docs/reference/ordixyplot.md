# Trellis (Lattice) Plots for Ordination

Function `ordixyplot` provides an interface to plot ordination results
using Trellis function
[`xyplot`](https://rdrr.io/pkg/lattice/man/xyplot.html) in package
lattice.

## Usage

``` r
ordixyplot(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi", aspect = "iso", envfit,
    type = c("p", "biplot"), ...)
```

## Arguments

- x:

  An ordination result that
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md)
  knows: any ordination result in vegan and many others.

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

- ...:

  Arguments passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md)
  methods or lattice functions.

## Details

The function provides an interface to the corresponding lattice
function. All graphical parameters are passed to the lattice function so
that these graphs are configurable. See
[`Lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html) and
[`xyplot`](https://rdrr.io/pkg/lattice/man/xyplot.html) for details,
usage and possibilities.

The argument `x` must always be an ordination result. The scores are
extracted with vegan function
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md) so
that these functions work with all vegan ordinations and many others.

The `formula` is used to define the models. Function has a simple
default formula which is used if `formula` is missing. The formula must
use the names of ordination scores and names of `data`.

The ordination scores are found from `x`, and `data` is optional. The
`data` should contain other variables than ordination scores to be used
in plots. Typically, they are environmental variables (typically
factors) to define panels or plot symbols.

The proper work is done by the panel function. The layout can be changed
by defining own panel functions. See
[`panel.xyplot`](https://rdrr.io/pkg/lattice/man/panel.xyplot.html) for
details and survey of possibilities.

Ordination graphics should always be isometric: same scale should be
used in all axes. This is controlled (and can be changed) with argument
`aspect` in `ordixyplot`.

## Value

The function return
[`Lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html) objects of
class `"trellis"`.

## Author

Jari Oksanen

## See also

[`Lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html),
[`xyplot`](https://rdrr.io/pkg/lattice/man/xyplot.html).

## Note

vegan releases 2.6-10 and earlier had lattice functions
[`ordicloud`](https://vegandevs.github.io/vegan/reference/vegan-defunct.md)
and
[`ordisplom`](https://vegandevs.github.io/vegan/reference/vegan-defunct.md)
which are now deprecated. However,
[vegan3d](https://CRAN.R-project.org/package=vegan3d) (version 1.4-0 and
later) has function `ordilattice3d` which is equal to `ordicloud`.

## Examples

``` r
data(dune, dune.env)
ord <- cca(dune)
## Scatter plot with polygons
ordixyplot(ord, data=dune.env, form = CA1 ~ CA2 | Management,
  groups=Manure, type = c("p","polygon"))

## Choose a different scaling
ordixyplot(ord, scaling = "sites")

## ... Slices of third axis
ordixyplot(ord, form = CA1 ~ CA2 | lattice::equal.count(CA3, 4),
   type = c("g","p", "polygon"))

## Display environmental variables
ordixyplot(ord, envfit = envfit(ord ~ Management + A1, dune.env, choices=1:3))
```
