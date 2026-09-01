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

`plot` methods for`renyiaccum` are provided by the
[ggvegan](https://CRAN.R-project.org/package=ggvegan) (`fortify`
function with an Example to design a graph).
[vegan3d](https://CRAN.R-project.org/package=vegan3d) can make dynamics
graphics with `rgl.renyiaccum`.

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
#> 14 4.584967 4.446502 4.302370 4.017494 3.570654 3.147989 2.879572 2.737063
#> 13 4.532599 4.402627 4.265383 3.982373 3.480484 2.941455 2.635488 2.495964
#> 40 4.382027 4.139704 3.855058 3.234849 2.450078 2.041264 1.838633 1.725318
#> 46 4.454347 4.298951 4.134690 3.810489 3.343102 2.975446 2.756061 2.629116
#> 23 4.595120 4.471984 4.339463 4.062575 3.588236 3.127012 2.864470 2.729487
#> 6  4.442651 4.279010 4.108089 3.776575 3.290256 2.843421 2.545504 2.385696
#> 32 4.477337 4.321957 4.152199 3.784873 3.135608 2.569002 2.295428 2.166841
#> 38 4.406719 4.210287 3.990829 3.516082 2.756942 2.203064 1.954323 1.834725
#> 20 4.605170 4.472551 4.336913 4.077327 3.683250 3.264270 2.924281 2.735886
#> 19 4.691348 4.542366 4.376816 4.013094 3.369175 2.810963 2.536968 2.413066
#> 24 4.553877 4.417824 4.273661 3.979427 3.487632 2.987362 2.659631 2.489169
#> 2  4.430817 4.280081 4.128601 3.848471 3.452678 3.084868 2.820363 2.656886
#>          32       64      Inf
#> 14 2.665029 2.626467 2.585711
#> 13 2.424969 2.387442 2.350154
#> 40 1.669834 1.643329 1.617652
#> 46 2.558311 2.519575 2.480266
#> 23 2.652997 2.611650 2.570849
#> 6  2.308989 2.272339 2.236834
#> 32 2.100247 2.067002 2.034706
#> 38 1.776031 1.747842 1.720532
#> 20 2.647690 2.605663 2.564949
#> 19 2.353743 2.322054 2.286548
#> 24 2.408956 2.370718 2.333676
#> 2  2.572204 2.531380 2.491827
mod <- renyiaccum(BCI[i,])
persp(mod)
```
