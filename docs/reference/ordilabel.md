# Add Text on Non-transparent Label to an Ordination Plot.

Function `ordilabel` is similar to
[`text`](https://rdrr.io/r/graphics/text.html), but the text is on an
opaque label. This can help in crowded ordination plots: you still
cannot see all text labels, but at least the uppermost ones are
readable. Argument `priority` helps to make the most important labels
visible. Function can be used in pipe after ordination `plot` or
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
command.

## Usage

``` r
ordilabel(x, display, labels, choices = c(1, 2), priority, select,
    cex = 0.8, fill = "white", border = NULL, col = NULL, xpd = TRUE, ...)
```

## Arguments

- x:

  An ordination object an any object known to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md).

- display:

  Kind of scores displayed (passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md)).

- labels:

  Optional text used in plots instead of the default. If `select` is
  given, the labels are given only to selected items in the order they
  occur in the scores.

- choices:

  Axes shown (passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md)).

- priority:

  Vector of the same length as the number of scores or selected items.
  The items with high priority will be plotted uppermost.

- select:

  Items to be displayed. This can either be a logical vector which is
  `TRUE` for displayed items or a vector of indices of displayed items.

- cex:

  Character expansion for the text (passed to
  [`text`](https://rdrr.io/r/graphics/text.html)).

- fill:

  Background colour of the labels (the `col` argument of
  [`polygon`](https://rdrr.io/r/graphics/polygon.html)).

- border:

  The colour and visibility of the border of the label as defined in
  [`polygon`](https://rdrr.io/r/graphics/polygon.html). The default is
  to use text colour `col`.

- col:

  Text colour.

- xpd:

  Draw labels also outside the plot region.

- ...:

  Other arguments (passed to
  [`text`](https://rdrr.io/r/graphics/text.html)).

## Details

The function may be useful with crowded ordination plots, in particular
together with argument `priority`. You will not see all text labels, but
at least some are readable. Function can be used as a part of a pipe
(`|>`) in place of `text` after an ordination `plot` command (see
Examples).

Other alternatives for cluttered plots are
[`identify.ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md),
[`orditorp`](https://vegandevs.github.io/vegan/reference/orditorp.md),
[`ordipointlabel`](https://vegandevs.github.io/vegan/reference/ordipointlabel.md),
and `orditkplot` ([vegan3d](https://CRAN.R-project.org/package=vegan3d)
package).

## Author

Jari Oksanen

## See also

[`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
and
[`text.ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
that can use the function with argument `bg`.

## Examples

``` r
data(dune)
ord <- cca(dune)
plot(ord, type = "n")
## add text
ordilabel(ord, dis="sites", cex=1.2, font=3, fill="hotpink", col="blue")
## You may prefer separate plots, but here species as well
ordilabel(ord, dis="sp", font=2, priority=colSums(dune))

## use in a pipe
plot(ord, type = "n") |>
   ordilabel("spec", font = 3, priority = colSums(dune)) |>
   points("sites", pch=21, bg = "yellow", col = "blue")
```
