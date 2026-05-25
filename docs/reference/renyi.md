# Renyi and Hill Diversities and Corresponding Accumulation Curves

Function `renyi` find Rényi diversities with any scale or the
corresponding Hill number (Hill 1973). Function `renyiaccum` finds these
statistics with accumulating sites.

## Usage

``` r
renyi(x, scales = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf),
   hill = FALSE)
renyiaccum(x, scales = c(0, 0.5, 1, 2, 4, Inf), permutations = 100,
    raw = FALSE, collector = FALSE, subset, ...)
# S3 method for class 'renyi'
plot(x, ...)
# S3 method for class 'renyiaccum'
plot(x, ...)
# S3 method for class 'renyiaccum'
lines(x, what, ...)
# S3 method for class 'renyiaccum'
persp(x, theta = 220, col = heat.colors(100), zlim, ...)
```

## Arguments

- x:

  Community data matrix or plotting object.

- scales:

  Scales of Rényi diversity.

- hill:

  Calculate Hill numbers.

- permutations:

  Usually an integer giving the number permutations, but can also be a
  list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or a
  permutation matrix where each row gives the permuted indices.

- raw:

  if `FALSE` then return summary statistics of permutations, and if
  `TRUE` then returns the individual permutations.

- collector:

  Accumulate the diversities in the order the sites are in the data set,
  and the collector curve can be plotted against summary of
  permutations. The argument is ignored if `raw = TRUE`.

- subset:

  logical expression indicating sites (rows) to keep: missing values are
  taken as `FALSE`.

- theta:

  Angle defining the viewing direction (azimuthal) in
  [`persp`](https://rdrr.io/r/graphics/persp.html).

- what:

  Lines to be drawn. These can be any item of dimension 3 of
  `renyiaccum` result, with defaults `"mean"`, `"min"`, `"max"`,
  `"Qnt 0.025"`, `"Qnt 0.975"`, and optional `"Collector"` or quoted
  numbers of permutations.

- col:

  Colours used for surface. Single colour will be passed on, and vector
  colours will be selected by the midpoint of a rectangle in
  [`persp`](https://rdrr.io/r/graphics/persp.html).

- zlim:

  Limits of vertical axis.

- ...:

  Other arguments which are passed to `renyi` and to graphical
  functions.

## Details

Common
[`diversity`](https://vegandevs.github.io/vegan/reference/diversity.md)
indices are special cases of Rényi diversity \$\$H_a = \frac{1}{1-a}
\log \sum p_i^a\$\$ where \\a\\ is a scale parameter, and Hill (1975)
suggested to use so-called ‘Hill numbers’ defined as \\N_a =
\exp(H_a)\\. Some Hill numbers are the number of species with \\a = 0\\,
\\\exp(H')\\ or the exponent of Shannon diversity with \\a = 1\\,
inverse Simpson with \\a = 2\\ and \\1/ \max(p_i)\\ with \\a = \infty\\.
According to the theory of diversity ordering, one community can be
regarded as more diverse than another only if its Rényi diversities are
all higher (Tóthmérész 1995).

Function `renyiaccum` is similar to
[`specaccum`](https://vegandevs.github.io/vegan/reference/specaccum.md)
but finds Rényi or Hill diversities at given `scales` for random
permutations of accumulated sites. It has a `persp` method to plot the
diversity surface against scale and number and sites.

`plot` and `lines` methods for `renyi` and `renyiaccum` are based on
[`matplot`](https://rdrr.io/r/graphics/matplot.html) and
[`matlines`](https://rdrr.io/r/graphics/matplot.html) and draw
everything on a single plot frame. In vegan pre version 2.8-0 these
graphs were based on lattice graphics that used separate panels. This
types of plots can be made in the
[ggvegan](https://CRAN.R-project.org/package=ggvegan) package (`fortify`
functions). [vegan3d](https://CRAN.R-project.org/package=vegan3d) can
make dynamics graphics with `rgl.renyiaccum`.

## Value

Function `renyi` returns a data frame of selected indices. Function
`renyiaccum` with argument `raw = FALSE` returns a three-dimensional
array, where the first dimension are the accumulated sites, second
dimension are the diversity scales, and third dimension are the summary
statistics `mean`, `stdev`, `min`, `max`, `Qnt 0.025` and `Qnt 0.975`.
With argument `raw = TRUE` the statistics on the third dimension are
replaced with individual permutation results.

## References

Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
consequences. *Ecology* 54, 427–473.

Kindt, R., Van Damme, P., Simons, A.J. (2006). Tree diversity in western
Kenya: using profiles to characterise richness and evenness.
*Biodiversity and Conservation* 15, 1253–1270.

Tóthmérész, B. (1995). Comparison of different methods for diversity
ordering. *Journal of Vegetation Science* 6, 283–290.

## Author

Roeland Kindt and Jari Oksanen

## See also

[`diversity`](https://vegandevs.github.io/vegan/reference/diversity.md)
for diversity indices, and
[`specaccum`](https://vegandevs.github.io/vegan/reference/specaccum.md)
for ordinary species accumulation curves. For alternative graphics, see
`fortify` function in
[ggvegan](https://CRAN.R-project.org/package=ggvegan) package.

## Examples

``` r
data(BCI)
i <- sample(nrow(BCI), 12)
mod <- renyi(BCI[i,])
plot(mod)

mod <- renyiaccum(BCI)
persp(mod)

## plot with empirical 95% confidence limits
plot(mod, lwd = 3)
lines(mod, "Qnt 0.025", lty=1)
lines(mod, "Qnt 0.975", lty=1)
```
