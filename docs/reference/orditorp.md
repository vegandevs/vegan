# Add Text or Points to Ordination Plots

The function adds [`text`](https://rdrr.io/r/graphics/text.html) or
[`points`](https://rdrr.io/r/graphics/points.html) to ordination plots.
Text will be used if this can be done without overwriting other text
labels, and points will be used otherwise. The function can help in
reducing clutter in ordination graphics, but manual editing may still be
necessary.

## Usage

``` r
orditorp(x, display, labels, choices = c(1, 2), priority,
    select, cex = 0.7, pcex, col = par("col"), pcol,
    pch = par("pch"), air = 1, ...)
```

## Arguments

- x:

  A result object from ordination or an
  [`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
  result. If the function is used in
  [`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
  pipe, this should be missing and first argument be `display`.

- display:

  Items to be displayed in the plot. Only one alternative is allowed.
  Typically this is `"sites"` or `"species"`.

- labels:

  Optional text used for labels. Row names of scores will be used if
  this is missing. If `select` is used, labels are given only selected
  items in the order they occur in the scores.

- choices:

  Axes shown.

- priority:

  Text will be used for items with higher priority if labels overlap.
  This should be vector of the same length as the number of items
  plotted or number of scores.

- select:

  Items to be displayed. This can either be a logical vector which is
  `TRUE` for displayed items or a vector of indices, or `labels` and if
  these are missing, row names of scores. If a logical vector is used,
  it must have the same length as the scores plotted.

- cex, pcex:

  Text and point sizes, see
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html)..

- col, pcol:

  Text and point colours, see
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html).

- pch:

  Plotting character, see
  [`points`](https://rdrr.io/r/graphics/points.html).

- air:

  Amount of empty space between text labels. Values \<1 allow
  overlapping text.

- ...:

  Other arguments to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md) (and
  its various methods), [`text`](https://rdrr.io/r/graphics/text.html)
  and [`points`](https://rdrr.io/r/graphics/points.html).

## Details

Function `orditorp` will add either text or points to an existing plot.
The items with high `priority` will be added first and
[`text`](https://rdrr.io/r/graphics/text.html) will be used if this can
be done without overwriting previous labels,and
[`points`](https://rdrr.io/r/graphics/points.html) will be used
otherwise. If `priority` is missing, labels will be added from the
outskirts to the centre. Function `orditorp` can be used with most
ordination results, or plotting results from
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md) or
ordination plot functions
([`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`plot.decorana`](https://vegandevs.github.io/vegan/reference/decorana.md),
[`plot.metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md)).
Function can also be used in a pipe (`|>`) where the first command is a
vegan ordination `plot` command or
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md).

Arguments can be passed to the relevant
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md) method
for the ordination object (`x`) being drawn. See the relevant
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md) help
page for arguments that can be used.

## Value

The function returns invisibly the The function returns invisibly a
logical vector where `TRUE` means that item was labelled with text and
`FALSE` means that it was marked with a point. If function is used in an
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
pipe, it will return the input `ordiplot` object, but amend the plotted
scores with this vector as attribute `"orditorp"`. The returned vector
can be used as the `select` argument in ordination `text` and `points`
functions.

## Author

Jari Oksanen

## Examples

``` r
## A cluttered ordination plot :
data(BCI)
mod <- cca(BCI)
plot(mod, dis="sp", type="t")

# Now with orditorp and abbreviated species names
cnam <- make.cepnames(names(BCI))
plot(mod, dis="sp", type="n")
stems <- colSums(BCI)
orditorp(mod, "sp", labels = cnam, priority=stems, pch="+", pcol="grey")


## show select in action
set.seed(1)
take <- sample(ncol(BCI), 50)
plot(mod, dis="sp", type="n")
stems <- colSums(BCI)
## only selected items are labelled, and the labels must be in the some
## order as in the scores
orditorp(mod, "sp", labels = cnam[sort(take)], priority=stems, select = take,
         pch="+", pcol="grey")

#> Warning: object '.Random.seed' not found
```
