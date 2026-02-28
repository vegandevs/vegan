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
#> 11  86 39.29329 19.38852 10.403636 6.082374 3.859814 2.636380 1.917208 1.467852
#> 15  92 41.60244 20.33533 10.818968 6.277031 3.956635 2.686972 1.944677 1.483199
#> 45  80 35.98439 17.64135  9.493581 5.607534 3.609518 2.502487 1.844372 1.427545
#> 24  94 42.55654 20.75558 10.992840 6.343926 3.979427 2.692702 1.944478 1.481404
#> 38  81 36.25129 17.62089  9.387396 5.495029 3.516082 2.432141 1.793969 1.392413
#> 3   89 39.79074 19.34183 10.294475 6.003596 3.814060 2.612130 1.905004 1.461939
#> 10  93 41.50369 20.10059 10.639881 6.164157 3.889803 2.648176 1.922207 1.470121
#> 44  80 36.52805 18.09206  9.782209 5.776652 3.705016 2.555595 1.873754 1.443791
#> 8   87 39.98200 19.77567 10.603913 6.181934 3.908381 2.659907 1.928610 1.473408
#> 47 101 44.24038 21.02130 10.943397 6.262025 3.920918 2.658316 1.925965 1.471937
#> 46  85 38.71044 19.06049 10.224553 5.987209 3.810489 2.611431 1.904913 1.461976
#> 35  82 33.26640 14.81511  7.386325 4.164622 2.641859 1.851529 1.400886 1.120527
#>         1.8         2
#> 11 1.171254 0.9658398
#> 15 1.180006 0.9709057
#> 45 1.148577 0.9528853
#> 24 1.178152 0.9694268
#> 38 1.124480 0.9365144
#> 3  1.168492 0.9646078
#> 10 1.172341 0.9663808
#> 44 1.157572 0.9578733
#> 8  1.173986 0.9671998
#> 47 1.173489 0.9672083
#> 46 1.168556 0.9646728
#> 35 0.932495 0.7983976
diversity(BCI[i,],"simpson") == x1[["2"]]
#>   11   15   45   24   38    3   10   44    8   47   46   35 
#> TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
x2 <- tsallis(BCI[i,],norm=TRUE)
x2
#>    0       0.2       0.4       0.6       0.8         1       1.2       1.4
#> 11 1 0.9081594 0.8567331 0.8377008 0.8430823 0.8642843 0.8927066 0.9212566
#> 15 1 0.9102071 0.8608373 0.8437404 0.8507238 0.8729284 0.9015588 0.9295308
#> 45 1 0.8821127 0.8163126 0.7912066 0.7963976 0.8213811 0.8559075 0.8914611
#> 24 1 0.9149445 0.8667059 0.8486326 0.8536824 0.8738548 0.9008891 0.9279024
#> 38 1 0.8797117 0.8089277 0.7777332 0.7771539 0.7978912 0.8304023 0.8662156
#> 3  1 0.8943582 0.8362208 0.8155352 0.8226868 0.8476046 0.8803722 0.9129172
#> 10 1 0.9001004 0.8450759 0.8255446 0.8324355 0.8561634 0.8872568 0.9180261
#> 44 1 0.8954397 0.8371683 0.8152612 0.8204162 0.8431126 0.8740718 0.9056624
#> 8  1 0.9154260 0.8674322 0.8491531 0.8535754 0.8729254 0.8992501 0.9258854
#> 47 1 0.8972246 0.8387068 0.8167095 0.8229465 0.8477709 0.8810155 0.9141216
#> 46 1 0.9032439 0.8485357 0.8278693 0.8331492 0.8554539 0.8856795 0.9162035
#> 35 1 0.7992527 0.6748192 0.6083778 0.5865659 0.5978626 0.6310844 0.6757387
#>          1.6       1.8         2
#> 11 0.9455721 0.9640735 0.9770705
#> 15 0.9527060 0.9698220 0.9814590
#> 45 0.9225835 0.9470176 0.9647964
#> 24 0.9507028 0.9678509 0.9797399
#> 38 0.8993700 0.9268720 0.9480763
#> 3  0.9403687 0.9610574 0.9754460
#> 10 0.9438798 0.9632975 0.9767719
#> 44 0.9330830 0.9544340 0.9698467
#> 8  0.9486740 0.9660687 0.9783170
#> 47 0.9418898 0.9625904 0.9767846
#> 46 0.9422698 0.9621112 0.9760219
#> 35 0.7233557 0.7684000 0.8081341
mod1 <- tsallisaccum(BCI[i,])
persp(mod1)

mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
```
