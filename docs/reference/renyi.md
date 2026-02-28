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

`plot` methods for `renyi` and `renyiaccum` are provided by the
[ggvegan](https://CRAN.R-project.org/package=ggvegan) (`autoplot`
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
for ordinary species accumulation curves.

## Examples

``` r
data(BCI)
i <- sample(nrow(BCI), 12)
renyi(BCI[i,])
#>           0     0.25      0.5        1        2        4        8       16
#> 4  4.543295 4.399850 4.253966 3.976563 3.561778 3.157132 2.878594 2.727145
#> 46 4.454347 4.298951 4.134690 3.810489 3.343102 2.975446 2.756061 2.629116
#> 43 4.454347 4.285742 4.104905 3.736254 3.145700 2.593930 2.282515 2.134836
#> 29 4.454347 4.290364 4.104978 3.688721 2.958515 2.387306 2.140083 2.028516
#> 8  4.477337 4.342394 4.199865 3.908381 3.417320 2.931933 2.644037 2.501526
#> 20 4.605170 4.472551 4.336913 4.077327 3.683250 3.264270 2.924281 2.735886
#> 6  4.442651 4.279010 4.108089 3.776575 3.290256 2.843421 2.545504 2.385696
#> 5  4.615121 4.460220 4.297238 3.969940 3.436618 2.885713 2.532071 2.365033
#> 39 4.430817 4.240081 4.021062 3.530494 2.749191 2.192896 1.918817 1.792300
#> 40 4.382027 4.139704 3.855058 3.234849 2.450078 2.041264 1.838633 1.725318
#> 33 4.454347 4.283915 4.102866 3.740392 3.186761 2.696061 2.399570 2.247627
#> 30 4.574711 4.408852 4.228938 3.851598 3.225546 2.661005 2.342443 2.189252
#>          32       64      Inf
#> 4  2.648728 2.607649 2.566920
#> 46 2.558311 2.519575 2.480266
#> 43 2.066044 2.033250 2.001480
#> 29 1.969346 1.938456 1.908170
#> 8  2.427926 2.389874 2.352536
#> 20 2.647690 2.605663 2.564949
#> 6  2.308989 2.272339 2.236834
#> 5  2.288749 2.252419 2.217225
#> 39 1.734488 1.706957 1.680286
#> 40 1.669834 1.643329 1.617652
#> 33 2.175329 2.140800 2.107350
#> 30 2.118649 2.085019 2.052441
mod <- renyiaccum(BCI[i,])
persp(mod)
```
