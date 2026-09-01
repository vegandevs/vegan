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
#> Run 1 stress 0.1955836 
#> Run 2 stress 0.1843196 
#> ... New best solution
#> ... Procrustes: rmse 7.106978e-06  max resid 2.214709e-05 
#> ... Similar to previous best
#> Run 3 stress 0.236479 
#> Run 4 stress 0.233173 
#> Run 5 stress 0.2380741 
#> Run 6 stress 0.2028828 
#> Run 7 stress 0.2079057 
#> Run 8 stress 0.1985583 
#> Run 9 stress 0.1845801 
#> ... Procrustes: rmse 0.04933527  max resid 0.1574168 
#> Run 10 stress 0.2048307 
#> Run 11 stress 0.1974406 
#> Run 12 stress 0.18584 
#> Run 13 stress 0.1955836 
#> Run 14 stress 0.2005511 
#> Run 15 stress 0.2097525 
#> Run 16 stress 0.2066172 
#> Run 17 stress 0.18458 
#> ... Procrustes: rmse 0.04934437  max resid 0.1574665 
#> Run 18 stress 0.2028828 
#> Run 19 stress 0.207032 
#> Run 20 stress 0.1825658 
#> ... New best solution
#> ... Procrustes: rmse 0.04161708  max resid 0.1517602 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     18: stress ratio > sratmax
#>      2: scale factor of the gradient < sfgrmin
stressplot(mod)

gof <- goodness(mod)
gof
#>  [1] 0.02984493 0.03513717 0.04189403 0.04598138 0.04003132 0.03441393
#>  [7] 0.03294968 0.03050088 0.03060814 0.02994098 0.03526199 0.02621444
#> [13] 0.03831043 0.02980902 0.03369408 0.02225920 0.03561632 0.03505253
#> [19] 0.06577463 0.03268362 0.03503005 0.02956645 0.05168185 0.04601906
plot(mod, display = "sites", type = "n")
points(mod, display = "sites", cex = 2*gof/mean(gof))
```
