# PCA biplot

Draws a PCA biplot with species scores indicated by biplot arrows

## Usage

``` r
# S3 method for class 'rda'
biplot(x, choices = c(1, 2), scaling = "species",
       display = c("sites", "species"), type, xlim, ylim, col = c(1,2), 
       const, correlation = FALSE, ...)
```

## Arguments

- x:

  A [`rda`](https://vegandevs.github.io/vegan/reference/cca.md) result
  object.

- choices:

  Axes to show.

- scaling:

  Scaling for species and site scores. Either species (`2`) or site
  (`1`) scores are scaled by eigenvalues, and the other set of scores is
  left unscaled, or with `3` both are scaled symmetrically by square
  root of eigenvalues. With negative scaling values in `rda`, species
  scores are divided by standard deviation of each species and
  multiplied with an equalizing constant. Unscaled raw scores stored in
  the result can be accessed with `scaling = 0`.

  The type of scores can also be specified as one of `"none"`,
  `"sites"`, `"species"`, or `"symmetric"`, which correspond to the
  values `0`, `1`, `2`, and `3` respectively. Argument `correlation` can
  be used in combination with these character descriptions to get the
  corresponding negative value.

- correlation:

  logical; if `scaling` is a character description of the scaling type,
  `correlation` can be used to select correlation-like scores for PCA.
  See argument `scaling` for details.

- display:

  Scores shown. These must some of the alternatives `"species"` for
  species scores, and/or `"sites"` for site scores.

- type:

  Type of plot: partial match to `text` for text labels, `points` for
  points, and `none` for setting frames only. If omitted, `text` is
  selected for smaller data sets, and `points` for larger. Can be of
  length 2 (e.g. `type = c("text", "points")`), in which case the first
  element describes how species scores are handled, and the second how
  site scores are drawn.

- xlim, ylim:

  the x and y limits (min, max) of the plot.

- col:

  Colours used for sites and species (in this order). If only one colour
  is given, it is used for both.

- const:

  General scaling constant for
  [`scores.rda`](https://vegandevs.github.io/vegan/reference/plot.cca.md).

- ...:

  Other parameters for plotting functions.

## Details

Produces a plot or biplot of the results of a call to
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md). It is
common for the "species" scores in a PCA to be drawn as biplot arrows
that point in the direction of increasing values for that variable. The
`biplot.rda` function provides a wrapper to `plot.cca` to allow the easy
production of such a plot.

`biplot.rda` is only suitable for unconstrained models. If used on an
ordination object with constraints, an error is issued.

Arrow heads are at the actual species scores, and the possible text is
projected to direction of the arrow.

## Value

The `plot` function returns invisibly a plotting structure which can be
used by
[`identify.ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
to identify the points or other functions in the
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
family.

## Author

Gavin Simpson and Jari Oksanen.

## Note

Prior to vegan 2.7-2 the species `text` was at the actual location and
arrows were shorter than with `points`.

## See also

[`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
can also draw biplot arrows since vegan 2.7-0.

## Examples

``` r
data(dune)
mod <- rda(dune, scale = TRUE)
biplot(mod, scaling = "symmetric")


## plot.cca can do the same
plot(mod, scaling = "symmetric", spe.par = list(arrows=TRUE))


## different type for species and site scores
biplot(mod, scaling = "symmetric", type = c("text", "points"))


## We can use ordiplot pipes to build similar plots with flexible
## control
plot(mod, scaling = "symmetric", type="n") |>
   points("sites", cex=0.7) |>
   text("species", arrows=TRUE, length=0.05, col=2, cex=0.7, font=3)
```
