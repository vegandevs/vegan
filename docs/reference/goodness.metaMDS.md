# Goodness of Fit and Shepard Plot for Nonmetric Multidimensional Scaling

Function `goodness.metaMDS` find goodness of fit measure for points in
nonmetric multidimensional scaling, and function `stressplot` makes a
[`Shepard`](https://rdrr.io/pkg/MASS/man/isoMDS.html) diagram.

## Usage

``` r
# S3 method for class 'metaMDS'
goodness(object, dis, ...)
# Default S3 method
stressplot(object, dis, pch, p.col = "blue", l.col = "red", 
    lwd = 2, ...)
```

## Arguments

- object:

  A result object from
  [`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md),
  [`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md) or
  [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html).

- dis:

  Dissimilarities. This should not be used with
  [`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md) or
  [`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md),
  but must be used with when the dissimilarities cannot be reconstructed
  from the result object.

- pch:

  Plotting character for points. Default is dependent on the number of
  points.

- p.col, l.col:

  Point and line colours.

- lwd:

  Line width. For
  [`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md)
  the default is `lwd = 1` if more than two lines are drawn, and
  `lwd = 2` otherwise.

- ...:

  Other parameters to functions, e.g. graphical parameters.

## Details

Function `goodness.metaMDS` finds a goodness of fit statistic for
observations (points). This is defined so that sum of squared values is
equal to squared stress. Large values indicate poor fit.

Function `stressplot` draws a Shepard diagram which is a plot of
ordination distances and monotone or linear fit line against original
dissimilarities. In addition, it displays two correlation-like
statistics on the goodness of fit in the graph. The nonmetric fit is
based on stress \\S\\ and defined as \\R^2 = 1-S^2\\. The “linear fit”
is the squared correlation between fitted values and ordination
distances. For
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md), the
“linear fit” and \\R^2\\ from “stress type 2” are equal.

Both functions can be used with
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md),
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md) and
[`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html). The original
dissimilarities should not be given for
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md) or
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md)
results, but they must given if the result object has no information to
reconstruct dissmilarities. The functions checks that dissimilarities
are consistent with current ordination, and refuses to analyse
inconsistent dissimilarities. Function `goodness.metaMDS` is generic in
vegan, but you must spell its name completely if the result has no
`class`.

## Value

Function `goodness` returns a vector of values. Function `stressplot`
returns invisibly an object with items for original dissimilarities,
ordination distances and fitted values.

## Author

Jari Oksanen.

## See also

[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md),
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md),
[`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html),
[`Shepard`](https://rdrr.io/pkg/MASS/man/isoMDS.html). Similar diagrams
for eigenvector ordinations can be drawn with
[`stressplot.wcmdscale`](https://vegandevs.github.io/vegan/reference/stressplot.wcmdscale.md),
[`stressplot.cca`](https://vegandevs.github.io/vegan/reference/stressplot.wcmdscale.md).

## Examples

``` r
data(varespec)
mod <- metaMDS(varespec)
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.1843196 
#> Run 1 stress 0.1967393 
#> Run 2 stress 0.2005512 
#> Run 3 stress 0.2178486 
#> Run 4 stress 0.2028828 
#> Run 5 stress 0.2414246 
#> Run 6 stress 0.2066172 
#> Run 7 stress 0.2390089 
#> Run 8 stress 0.1825658 
#> ... New best solution
#> ... Procrustes: rmse 0.04162267  max resid 0.1517847 
#> Run 9 stress 0.2403423 
#> Run 10 stress 0.2085949 
#> Run 11 stress 0.1825658 
#> ... New best solution
#> ... Procrustes: rmse 4.563076e-06  max resid 1.492218e-05 
#> ... Similar to previous best
#> Run 12 stress 0.2415059 
#> Run 13 stress 0.196245 
#> Run 14 stress 0.2265716 
#> Run 15 stress 0.222519 
#> Run 16 stress 0.2511508 
#> Run 17 stress 0.2097525 
#> Run 18 stress 0.18458 
#> Run 19 stress 0.2175652 
#> Run 20 stress 0.2048307 
#> *** Best solution repeated 1 times
stressplot(mod)

gof <- goodness(mod)
gof
#>  [1] 0.02984517 0.03513714 0.04189226 0.04598241 0.04003135 0.03441441
#>  [7] 0.03294885 0.03050048 0.03060771 0.02994096 0.03526275 0.02621428
#> [13] 0.03831060 0.02980905 0.03369533 0.02225870 0.03561596 0.03505236
#> [19] 0.06577470 0.03268424 0.03503075 0.02956664 0.05168036 0.04601986
plot(mod, display = "sites", type = "n")
points(mod, display = "sites", cex = 2*gof/mean(gof))
```
