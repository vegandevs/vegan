# Ordination Plots with Points and Optimized Locations for Text

Function produces ordination plots with labelled points. The points are
in the fixed locations given by the ordination, but the locations of the
text labels are optimized to minimize overplotting. The function is
useful with moderately crowded ordination plots.

## Usage

``` r
ordipointlabel(x, display = c("sites", "species"), choices = c(1, 2),
   col = c(1, 2),  pch = c("o", "+"), font = c(1, 1), 
   cex = c(0.7, 0.7), add = inherits(x, "ordiplot"), labels, bg, select, ...)

# S3 method for class 'ordipointlabel'
plot(x, ...)
```

## Arguments

- x:

  For `ordipointlabel` a result object from an ordination function or an
  ordination plot (possibly in a pipe). For `plot.ordipointlabel` an
  object from `ordipointlabel`.

- display:

  Scores displayed in the plot. The default is to show `"sites"` and
  `"species"` that are available to many ordination methods, but there
  can be only one set or more than two set of scores to `display`.

- choices:

  Axes shown.

- col, pch, font, cex:

  Colours, point types, font style and character expansion for each kind
  of scores displayed in the plot. These should be vectors of the same
  length as the number of items in `display`, or if there is only one
  `display` they can be a vector of the length of number items.

- add:

  Add to an existing plot. Default is `add = TRUE` when the function is
  used in a pipe, and `FALSE` usually.

- labels:

  Labels used in graph. Species (variable) and SU (row) names are used
  if this is missing. Labels must be given in one vector for all scores
  of `display`. Function `labels` can extract the current name from a
  saved `ordipointlabel` object. If `select` is used, labels are given
  only for the selected items.

- bg:

  Background colour for labels. If this is given, texts is drawn over
  non-transparent background. Either a single colour or vector of
  colours for each `display`, or with one display, for each label.

- select:

  Items to be displayed. This can either be a logical vector which is
  `TRUE` for displayed items or a vector of indices of displayed items.
  `select` is only used if a single set of scores is being plotted (i.e.
  `length(display) == 1`), otherwise it is ignored and a warning issued.
  If a logical vector is used, it must have the same length as the
  scores plotted.

- ...:

  Other arguments passed to
  [`points`](https://rdrr.io/r/graphics/points.html) and
  [`text`](https://rdrr.io/r/graphics/text.html).

## Details

The function uses simulated annealing
([`optim`](https://rdrr.io/r/stats/optim.html), `method = "SANN"`) to
optimize the locations of the text labels to the points. There are eight
possible locations: up, down, two sides and four corners. There is a
weak preference to text away from zero, and a weak avoidance of corners.
The locations and goodness of solution varies between runs, and there is
no guarantee of finding the global optimum, or the same text locations
twice. The optimization can take a long time in difficult cases with a
high number of potential overlaps. Several sets of scores can be
displayed in one plot.

The function can be used in a pipe where the first command is an
ordination `plot` command with `type = "n"` or to add points and lablels
to save vegan ordination plot object. See examples.

## Value

The function returns invisibly an object of class `ordipointlabel` with
items `xy` for coordinates of points, `labels` for coordinates of
labels, items `pch`, `cex` and `font` for graphical parameters of each
point or label. In addition, it returns the result of
[`optim`](https://rdrr.io/r/stats/optim.html) as an attribute `"optim"`.
The unit of overlap is the area of character `"m"`, and with varying
graphical parameters the smallest alternative.

There is a `plot` method based on `orditkplot` but it does not alter or
reset the graphical parameters via `par`.

The result object from `ordipointlabel` is similar as the `orditkplot`
object of [vegan3d](https://CRAN.R-project.org/package=vegan3d) package,
and it may be possible to further edit the result object with
`orditkplot`, but for good results it is necessary that the points span
the whole horizontal axis without empty margins.

## See also

The function is invoked for one set of scores (one `display`) from
[`text.ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
and
[`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
with argument `optimize = TRUE`.

## Author

Jari Oksanen

## Examples

``` r
data(dune, dune.env)
ord <- cca(dune)
ordipointlabel(ord)


## Use in a pipe: optimize species, sites & centroids together
ord <- cca(dune ~ Management + Moisture, dune.env)
plot(ord, scaling = "symmetric", type = "n") |>
  ordipointlabel(c("sites","species","centroids"), cex=c(0.7,0.7,1),
      col = c("black","red","blue"), font = c(1,3,1), pch=c(1,3,4), xpd=TRUE) |>
  text("biplot", col = "blue", bg = "white", cex=1)
```
