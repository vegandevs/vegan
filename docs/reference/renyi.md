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
#> 25 4.653960 4.516361 4.370520 4.074718 3.597768 3.151248 2.869456 2.703842
#> 2  4.430817 4.280081 4.128601 3.848471 3.452678 3.084868 2.820363 2.656886
#> 12 4.430817 4.264581 4.081595 3.698414 3.102424 2.609704 2.341802 2.209194
#> 7  4.406719 4.262593 4.115994 3.836811 3.417368 3.010718 2.720940 2.552329
#> 47 4.624973 4.458358 4.278871 3.920918 3.417580 3.043219 2.837600 2.715794
#> 13 4.532599 4.402627 4.265383 3.982373 3.480484 2.941455 2.635488 2.495964
#> 8  4.477337 4.342394 4.199865 3.908381 3.417320 2.931933 2.644037 2.501526
#> 35 4.418841 4.071199 3.633626 2.641859 1.601458 1.178736 1.016178 0.948487
#> 46 4.454347 4.298951 4.134690 3.810489 3.343102 2.975446 2.756061 2.629116
#> 17 4.532599 4.360865 4.164991 3.736897 3.090319 2.664944 2.463809 2.345745
#> 4  4.543295 4.399850 4.253966 3.976563 3.561778 3.157132 2.878594 2.727145
#> 1  4.532599 4.398890 4.265075 4.018412 3.674161 3.371787 3.167106 3.042194
#>           32       64       Inf
#> 25 2.6177623 2.576215 2.5359618
#> 2  2.5722044 2.531380 2.4918271
#> 12 2.1408236 2.106912 2.0739919
#> 7  2.4702953 2.431085 2.3930991
#> 47 2.6377435 2.596146 2.5555816
#> 13 2.4249694 2.387442 2.3501535
#> 8  2.4279255 2.389874 2.3525360
#> 35 0.9178907 0.903321 0.8892066
#> 46 2.5583114 2.519575 2.4802663
#> 17 2.2760275 2.240131 2.2051298
#> 4  2.6487282 2.607649 2.5669198
#> 1  2.9707634 2.930598 2.8859174
mod <- renyiaccum(BCI[i,])
persp(mod)
```
