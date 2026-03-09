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
#> 7   81 37.47840 18.73414 10.172688 6.004908 3.836811 2.631757 1.918091 1.469869
#> 5  100 44.20598 21.17686 11.081257 6.350040 3.969940 2.683412 1.937940 1.477217
#> 27  98 43.84660 21.16041 11.110930 6.371476 3.980281 2.687154 1.938584 1.476662
#> 26  90 40.93055 20.10590 10.740144 6.250173 3.947749 2.684219 1.943942 1.483076
#> 9   89 39.96271 19.39982 10.272025 5.952153 3.761331 2.568371 1.872010 1.438318
#> 46  85 38.71044 19.06049 10.224553 5.987209 3.810489 2.611431 1.904913 1.461976
#> 24  94 42.55654 20.75558 10.992840 6.343926 3.979427 2.692702 1.944478 1.481404
#> 2   83 38.06254 18.90957 10.230614 6.028061 3.848471 2.638667 1.922477 1.472692
#> 38  81 36.25129 17.62089  9.387396 5.495029 3.516082 2.432141 1.793969 1.392413
#> 25 104 46.16214 22.10667 11.522497 6.562459 4.074718 2.736660 1.965873 1.492328
#> 14  97 43.58438 21.14034 11.153293 6.418569 4.017494 2.713442 1.956269 1.488283
#>         1.8         2
#> 42 1.183858 0.9731849
#> 7  1.173083 0.9672014
#> 5  1.175553 0.9678267
#> 27 1.174682 0.9669962
#> 26 1.180039 0.9709567
#> 9  1.152087 0.9534257
#> 46 1.168556 0.9646728
#> 24 1.178152 0.9694268
#> 2  1.174887 0.9683393
#> 38 1.124480 0.9365144
#> 25 1.183958 0.9726152
#> 14 1.182226 0.9718626
diversity(BCI[i,],"simpson") == x1[["2"]]
#>   42    7    5   27   26    9   46   24    2   38   25   14 
#> TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
x2 <- tsallis(BCI[i,],norm=TRUE)
x2
#>    0       0.2       0.4       0.6       0.8         1       1.2       1.4
#> 42 1 0.9185764 0.8746393 0.8602787 0.8675785 0.8881987 0.9141183 0.9390681
#> 7  1 0.9094899 0.8600339 0.8427936 0.8492653 0.8706729 0.8985570 0.9261475
#> 5  1 0.9038026 0.8502581 0.8308757 0.8372445 0.8602030 0.8904867 0.9204832
#> 27 1 0.9112896 0.8605499 0.8410587 0.8456725 0.8661975 0.8940909 0.9221774
#> 26 1 0.9116552 0.8630994 0.8463498 0.8532945 0.8751655 0.9033025 0.9307642
#> 9  1 0.8982235 0.8387276 0.8137568 0.8156375 0.8358867 0.8656242 0.8971059
#> 46 1 0.9032439 0.8485357 0.8278693 0.8331492 0.8554539 0.8856795 0.9162035
#> 24 1 0.9149445 0.8667059 0.8486326 0.8536824 0.8738548 0.9008891 0.9279024
#> 2  1 0.9055048 0.8546877 0.8377974 0.8455660 0.8685692 0.8978637 0.9264249
#> 38 1 0.8797117 0.8089277 0.7777332 0.7771539 0.7978912 0.8304023 0.8662156
#> 25 1 0.9142073 0.8658106 0.8481902 0.8541991 0.8755378 0.9035478 0.9310603
#> 14 1 0.9134203 0.8653450 0.8483542 0.8548076 0.8762317 0.9040571 0.9313075
#>          1.6       1.8         2
#> 42 0.9595063 0.9744477 0.9845010
#> 7  0.9493990 0.9669335 0.9791421
#> 5  0.9456409 0.9644767 0.9775050
#> 27 0.9460502 0.9641597 0.9768635
#> 26 0.9535106 0.9703141 0.9817451
#> 9  0.9251745 0.9475646 0.9641383
#> 46 0.9422698 0.9621112 0.9760219
#> 24 0.9507028 0.9678509 0.9797399
#> 2  0.9501791 0.9678597 0.9800060
#> 38 0.8993700 0.9268720 0.9480763
#> 25 0.9538440 0.9706141 0.9819673
#> 14 0.9538903 0.9705576 0.9818818
mod1 <- tsallisaccum(BCI[i,])
persp(mod1)

mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
```
