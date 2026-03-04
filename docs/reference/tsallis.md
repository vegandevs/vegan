# Tsallis Diversity and Corresponding Accumulation Curves

Function `tsallis` find Tsallis diversities with any scale or the
corresponding evenness measures. Function `tsallisaccum` finds these
statistics with accumulating sites.

## Usage

``` r
tsallis(x, scales = seq(0, 2, 0.2), norm = FALSE, hill = FALSE)
tsallisaccum(x, scales = seq(0, 2, 0.2), permutations = 100, 
   raw = FALSE, subset, ...)
# S3 method for class 'tsallisaccum'
persp(x, theta = 220, phi = 15, col = heat.colors(100), zlim, ...)
```

## Arguments

- x:

  Community data matrix or plotting object.

- scales:

  Scales of Tsallis diversity.

- norm:

  Logical, if `TRUE` diversity values are normalized by their maximum
  (diversity value at equiprobability conditions).

- hill:

  Calculate Hill numbers.

- permutations:

  Usually an integer giving the number permutations, but can also be a
  list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or a
  permutation matrix where each row gives the permuted indices.

- raw:

  If `FALSE` then return summary statistics of permutations, and if TRUE
  then returns the individual permutations.

- subset:

  logical expression indicating sites (rows) to keep: missing values are
  taken as `FALSE`.

- theta, phi:

  angles defining the viewing direction. `theta` gives the azimuthal
  direction and `phi` the colatitude.

- col:

  Colours used for surface.

- zlim:

  Limits of vertical axis.

- ...:

  Other arguments which are passed to `tsallis` and to graphical
  functions.

## Details

The Tsallis diversity (also equivalent to Patil and Taillie diversity)
is a one-parametric generalised entropy function, defined as:

\$\$H_q = \frac{1}{q-1} (1-\sum\_{i=1}^S p_i^q)\$\$

where \\q\\ is a scale parameter, \\S\\ the number of species in the
sample (Tsallis 1988, Tothmeresz 1995). This diversity is concave for
all \\q\>0\\, but non-additive (Keylock 2005). For \\q=0\\ it gives the
number of species minus one, as \\q\\ tends to 1 this gives Shannon
diversity, for \\q=2\\ this gives the Simpson index (see function
[`diversity`](https://vegandevs.github.io/vegan/reference/diversity.md)).

If `norm = TRUE`, `tsallis` gives values normalized by the maximum:

\$\$H_q(max) = \frac{S^{1-q}-1}{1-q}\$\$

where \\S\\ is the number of species. As \\q\\ tends to 1, maximum is
defined as \\ln(S)\\.

If `hill = TRUE`, `tsallis` gives Hill numbers (numbers equivalents, see
Jost 2007):

\$\$D_q = (1-(q-1) H)^{1/(1-q)}\$\$

Details on plotting methods and accumulating values can be found on the
help pages of the functions
[`renyi`](https://vegandevs.github.io/vegan/reference/renyi.md) and
[`renyiaccum`](https://vegandevs.github.io/vegan/reference/renyi.md).

## Value

Function `tsallis` returns a data frame of selected indices. Function
`tsallisaccum` with argument `raw = FALSE` returns a three-dimensional
array, where the first dimension are the accumulated sites, second
dimension are the diversity scales, and third dimension are the summary
statistics `mean`, `stdev`, `min`, `max`, `Qnt 0.025` and `Qnt 0.975`.
With argument `raw = TRUE` the statistics on the third dimension are
replaced with individual permutation results.

## References

Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs
statistics. *J. Stat. Phis.* 52, 479–487.

Tothmeresz, B. (1995) Comparison of different methods for diversity
ordering. *Journal of Vegetation Science* **6**, 283–290.

Patil, G. P. and Taillie, C. (1982) Diversity as a concept and its
measurement. *J. Am. Stat. Ass.* **77**, 548–567.

Keylock, C. J. (2005) Simpson diversity and the Shannon-Wiener index as
special cases of a generalized entropy. *Oikos* **109**, 203–207.

Jost, L (2007) Partitioning diversity into independent alpha and beta
components. *Ecology* **88**, 2427–2439.

## Author

Péter Sólymos, <solymos@ualberta.ca>, based on the code of Roeland Kindt
and Jari Oksanen written for `renyi`

## See also

Accumulation routines are based on functions
[`renyi`](https://vegandevs.github.io/vegan/reference/renyi.md) and
[`renyiaccum`](https://vegandevs.github.io/vegan/reference/renyi.md),
and plots are provided by
[ggvegan](https://CRAN.R-project.org/package=ggvegan). An object of
class `tsallisaccum` can be displayed with dynamic 3D function
`rgl.renyiaccum` in the vegan3d package. See also settings for
[`persp`](https://rdrr.io/r/graphics/persp.html).

## Examples

``` r
data(BCI)
i <- sample(nrow(BCI), 12)
x1 <- tsallis(BCI[i,])
x1
#>      0      0.2      0.4       0.6      0.8        1      1.2      1.4      1.6
#> 42  86 39.74400 19.79375 10.684037 6.259100 3.966614 2.699614 1.954275 1.489483
#> 2   83 38.06254 18.90957 10.230614 6.028061 3.848471 2.638667 1.922477 1.472692
#> 5  100 44.20598 21.17686 11.081257 6.350040 3.969940 2.683412 1.937940 1.477217
#> 22  90 40.16368 19.41416 10.256297 5.939677 3.755413 2.566788 1.872586 1.439722
#> 50  92 41.30206 20.08829 10.663128 6.187254 3.906616 2.659292 1.929248 1.474476
#> 3   89 39.79074 19.34183 10.294475 6.003596 3.814060 2.612130 1.905004 1.461939
#> 9   89 39.96271 19.39982 10.272025 5.952153 3.761331 2.568371 1.872010 1.438318
#> 11  86 39.29329 19.38852 10.403636 6.082374 3.859814 2.636380 1.917208 1.467852
#> 18  88 40.68467 20.15164 10.784931 6.264583 3.944985 2.675802 1.935407 1.476270
#> 44  80 36.52805 18.09206  9.782209 5.776652 3.705016 2.555595 1.873754 1.443791
#> 8   87 39.98200 19.77567 10.603913 6.181934 3.908381 2.659907 1.928610 1.473408
#> 19 108 47.27018 22.32602 11.498387 6.492802 4.013094 2.692466 1.936696 1.473828
#>         1.8         2
#> 42 1.183858 0.9731849
#> 2  1.174887 0.9683393
#> 5  1.175553 0.9678267
#> 22 1.153640 0.9548316
#> 50 1.174993 0.9679784
#> 3  1.168492 0.9646078
#> 9  1.152087 0.9534257
#> 11 1.171254 0.9658398
#> 18 1.175165 0.9676685
#> 44 1.157572 0.9578733
#> 8  1.173986 0.9671998
#> 19 1.172481 0.9655820
diversity(BCI[i,],"simpson") == x1[["2"]]
#>   42    2    5   22   50    3    9   11   18   44    8   19 
#> TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
x2 <- tsallis(BCI[i,],norm=TRUE)
x2
#>    0       0.2       0.4       0.6       0.8         1       1.2       1.4
#> 42 1 0.9185764 0.8746393 0.8602787 0.8675785 0.8881987 0.9141183 0.9390681
#> 2  1 0.9055048 0.8546877 0.8377974 0.8455660 0.8685692 0.8978637 0.9264249
#> 5  1 0.9038026 0.8502581 0.8308757 0.8372445 0.8602030 0.8904867 0.9204832
#> 22 1 0.8945746 0.8334045 0.8082215 0.8109045 0.8325271 0.8637842 0.8965990
#> 50 1 0.9036351 0.8503796 0.8315870 0.8385564 0.8618931 0.8922712 0.9221556
#> 3  1 0.8943582 0.8362208 0.8155352 0.8226868 0.8476046 0.8803722 0.9129172
#> 9  1 0.8982235 0.8387276 0.8137568 0.8156375 0.8358867 0.8656242 0.8971059
#> 11 1 0.9081594 0.8567331 0.8377008 0.8430823 0.8642843 0.8927066 0.9212566
#> 18 1 0.9228941 0.8775173 0.8589786 0.8616917 0.8788828 0.9032164 0.9283108
#> 44 1 0.8954397 0.8371683 0.8152612 0.8204162 0.8431126 0.8740718 0.9056624
#> 8  1 0.9154260 0.8674322 0.8491531 0.8535754 0.8729254 0.8992501 0.9258854
#> 19 1 0.9079031 0.8537683 0.8315799 0.8347890 0.8554245 0.8846675 0.9147428
#>          1.6       1.8         2
#> 42 0.9595063 0.9744477 0.9845010
#> 2  0.9501791 0.9678597 0.9800060
#> 5  0.9456409 0.9644767 0.9775050
#> 22 0.9256369 0.9486075 0.9654409
#> 50 0.9471031 0.9657020 0.9784999
#> 3  0.9403687 0.9610574 0.9754460
#> 9  0.9251745 0.9475646 0.9641383
#> 11 0.9455721 0.9640735 0.9770705
#> 18 0.9500471 0.9667897 0.9786648
#> 44 0.9330830 0.9544340 0.9698467
#> 8  0.9486740 0.9660687 0.9783170
#> 19 0.9406572 0.9605041 0.9745226
mod1 <- tsallisaccum(BCI[i,])
persp(mod1)

mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
```
