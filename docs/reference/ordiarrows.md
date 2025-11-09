# Add Arrows and Line Segments to Ordination Diagrams

Functions to add arrows, line segments, regular grids of points. The
ordination diagrams can be produced by `vegan`
[`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`plot.decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)
or
[`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md).

## Usage

``` r
ordiarrows(ord, groups, levels, replicates, order.by, display = "sites",
         col = 1, show.groups, startmark, label = FALSE, length = 0.1, ...)
ordisegments(ord, groups, levels, replicates, order.by, display = "sites",
         col = 1, show.groups, label = FALSE, ...)
ordigrid(ord, levels, replicates, display = "sites",  lty = c(1,1), 
         col = c(1,1), lwd = c(1,1), ...)
```

## Arguments

- ord:

  An ordination object or an
  [`ordiplot`](https://vegandevs.github.io/vegan/reference/ordiplot.md)
  object.

- groups:

  Factor giving the groups for which the graphical item is drawn.

- levels, replicates:

  Alternatively, regular groups can be defined with arguments `levels`
  and `replicates`, where `levels` gives the number of groups, and
  `replicates` the number of successive items at the same group.

- order.by:

  Order points by increasing order of this variable within `groups`.
  Reverse sign of the variable for decreasing ordering.

- display:

  Item to displayed.

- show.groups:

  Show only given groups. This can be a vector, or `TRUE` if you want to
  show items for which condition is `TRUE`. This argument makes it
  possible to use different colours and line types for groups. The
  default is to show all groups.

- label:

  Label the `groups` by their names. In `ordiellipse`, `ordihull` and
  `ordispider` the the group name is in the centroid of the object, in
  `ordiarrows` in the start of the arrow, and in `ordisegments` at both
  ends. `ordiellipse` and `ordihull` use standard
  [`text`](https://rdrr.io/r/graphics/text.html), and others use
  [`ordilabel`](https://vegandevs.github.io/vegan/reference/ordilabel.md).

- startmark:

  plotting character used to mark the first item. The default is to use
  no mark, and for instance, `startmark = 1` will draw a circle. For
  other plotting characters, see `pch` in
  [`points`](https://rdrr.io/r/graphics/points.html).

- col:

  Colour of lines, `label` borders and `startmark` in `ordiarrows` and
  `ordisegments`. This can be a vector recycled for `groups`. In
  `ordigrid` it can be a vector of length 2 used for `levels` and
  `replicates`.

- length:

  Length of edges of the arrow head (in inches).

- lty, lwd:

  Line type, line width used for `level`s and `replicate`s in
  `ordigrid`.

- ...:

  Parameters passed to graphical functions such as
  [`lines`](https://rdrr.io/r/graphics/lines.html),
  [`segments`](https://rdrr.io/r/graphics/segments.html),
  [`arrows`](https://rdrr.io/r/graphics/arrows.html), or to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md) to
  select axes and scaling etc.

## Details

Function `ordiarrows` draws
[`arrows`](https://rdrr.io/r/graphics/arrows.html) and `ordisegments`
draws line [`segments`](https://rdrr.io/r/graphics/segments.html)
between successive items in the groups. Function `ordigrid` draws line
[`segments`](https://rdrr.io/r/graphics/segments.html) both within the
groups and for the corresponding items among the groups.

## Note

These functions add graphical items to ordination graph: You must draw a
graph first.

## Author

Jari Oksanen

## See also

The functions pass parameters to basic graphical functions, and you may
wish to change the default values in
[`arrows`](https://rdrr.io/r/graphics/arrows.html),
[`lines`](https://rdrr.io/r/graphics/lines.html) and
[`segments`](https://rdrr.io/r/graphics/segments.html). You can pass
parameters to
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md) as
well.

## Examples

``` r
example(pyrifos)
#> 
#> pyrifs> data(pyrifos)
#> 
#> pyrifs> ditch <- gl(12, 1, length=132)
#> 
#> pyrifs> week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
#> 
#> pyrifs> dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
mod <- rda(pyrifos)
plot(mod, type = "n")
## Annual succession by ditches, colour by dose
ordiarrows(mod, ditch, label = TRUE, col = as.numeric(dose))
legend("topright", levels(dose), lty=1, col=1:5, title="Dose")

## Show only control and highest Pyrifos treatment
plot(mod, type = "n")
ordiarrows(mod, ditch, label = TRUE, 
   show.groups = c("2", "3", "5", "11"))
ordiarrows(mod, ditch, label = TRUE, show = c("6", "9"),
   col = 2)
legend("topright", c("Control", "Pyrifos 44"), lty = 1, col = c(1,2))
```
