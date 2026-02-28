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
#> 22 4.510860 4.339960 4.152636 3.755413 3.097358 2.535617 2.246595 2.106328
#> 19 4.691348 4.542366 4.376816 4.013094 3.369175 2.810963 2.536968 2.413066
#> 4  4.543295 4.399850 4.253966 3.976563 3.561778 3.157132 2.878594 2.727145
#> 34 4.521789 4.360121 4.185237 3.821669 3.238761 2.748477 2.486674 2.352647
#> 13 4.532599 4.402627 4.265383 3.982373 3.480484 2.941455 2.635488 2.495964
#> 39 4.430817 4.240081 4.021062 3.530494 2.749191 2.192896 1.918817 1.792300
#> 38 4.406719 4.210287 3.990829 3.516082 2.756942 2.203064 1.954323 1.834725
#> 45 4.394449 4.203234 4.001504 3.609518 3.055171 2.600948 2.332329 2.189116
#> 47 4.624973 4.458358 4.278871 3.920918 3.417580 3.043219 2.837600 2.715794
#> 42 4.465908 4.336864 4.207514 3.966614 3.618790 3.295911 3.065074 2.923061
#> 2  4.430817 4.280081 4.128601 3.848471 3.452678 3.084868 2.820363 2.656886
#> 29 4.454347 4.290364 4.104978 3.688721 2.958515 2.387306 2.140083 2.028516
#>          32       64      Inf
#> 22 2.038748 2.006388 1.975038
#> 19 2.353743 2.322054 2.286548
#> 4  2.648728 2.607649 2.566920
#> 34 2.281821 2.245829 2.210738
#> 13 2.424969 2.387442 2.350154
#> 39 1.734488 1.706957 1.680286
#> 38 1.776031 1.747842 1.720532
#> 45 2.118990 2.085357 2.052773
#> 47 2.637744 2.596146 2.555582
#> 42 2.848959 2.810462 2.767769
#> 2  2.572204 2.531380 2.491827
#> 29 1.969346 1.938456 1.908170
mod <- renyiaccum(BCI[i,])
persp(mod)
```
