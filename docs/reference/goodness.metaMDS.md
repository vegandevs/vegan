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
  but must be used with
  [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html).

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
equal to squared stress. Large values indicate poor fit. The absolute
values of the goodness statistic depend on the definition of the stress:
[`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html) expresses stress in
percents, and therefore its goodness values are 100 times higher than
those of
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md)
which expresses the stress as a proportion.

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
results (the latter tries to reconstruct the dissimilarities using
[`metaMDSredist`](https://vegandevs.github.io/vegan/reference/metaMDS.md)
if [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html) was used as its
engine). With [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html) the
dissimilarities must be given. In either case, the functions inspect
that dissimilarities are consistent with current ordination, and refuse
to analyse inconsistent dissimilarities. Function `goodness.metaMDS` is
generic in vegan, but you must spell its name completely with
[`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html) which has no class.

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
#> Run 1 stress 0.2093084 
#> Run 2 stress 0.1985581 
#> Run 3 stress 0.2069724 
#> Run 4 stress 0.1825658 
#> ... New best solution
#> ... Procrustes: rmse 0.04162213  max resid 0.15177 
#> Run 5 stress 0.1825658 
#> ... New best solution
#> ... Procrustes: rmse 1.544563e-05  max resid 5.216924e-05 
#> ... Similar to previous best
#> Run 6 stress 0.195049 
#> Run 7 stress 0.18584 
#> Run 8 stress 0.2088293 
#> Run 9 stress 0.2218502 
#> Run 10 stress 0.2239801 
#> Run 11 stress 0.1825658 
#> ... Procrustes: rmse 2.359921e-05  max resid 6.590669e-05 
#> ... Similar to previous best
#> Run 12 stress 0.2079059 
#> Run 13 stress 0.18458 
#> Run 14 stress 0.1843196 
#> Run 15 stress 0.2287525 
#> Run 16 stress 0.1967393 
#> Run 17 stress 0.2069724 
#> Run 18 stress 0.1825658 
#> ... Procrustes: rmse 4.755455e-06  max resid 1.256699e-05 
#> ... Similar to previous best
#> Run 19 stress 0.2335413 
#> Run 20 stress 0.2388082 
#> *** Best solution repeated 3 times
stressplot(mod)

gof <- goodness(mod)
gof
#>  [1] 0.02984518 0.03513708 0.04189069 0.04598289 0.04003135 0.03441482
#>  [7] 0.03294830 0.03050074 0.03060725 0.02994076 0.03526393 0.02621416
#> [13] 0.03831048 0.02980900 0.03369686 0.02225850 0.03561549 0.03505228
#> [19] 0.06577482 0.03268427 0.03503154 0.02956646 0.05167921 0.04602049
plot(mod, display = "sites", type = "n")
points(mod, display = "sites", cex = 2*gof/mean(gof))
```
