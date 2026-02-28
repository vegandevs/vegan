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
#> 21  98 43.64355 21.01620 11.036849 6.340108 3.969925 2.685700 1.940242 1.478977
#> 35  82 33.26640 14.81511  7.386325 4.164622 2.641859 1.851529 1.400886 1.120527
#> 33  85 38.33236 18.75302 10.027610 5.869003 3.740392 2.569560 1.879580 1.446458
#> 5  100 44.20598 21.17686 11.081257 6.350040 3.969940 2.683412 1.937940 1.477217
#> 16  92 41.45601 20.19256 10.715725 6.210870 3.916821 2.663870 1.931567 1.475861
#> 7   81 37.47840 18.73414 10.172688 6.004908 3.836811 2.631757 1.918091 1.469869
#> 14  97 43.58438 21.14034 11.153293 6.418569 4.017494 2.713442 1.956269 1.488283
#> 22  90 40.16368 19.41416 10.256297 5.939677 3.755413 2.566788 1.872586 1.439722
#> 48  90 40.88736 20.03984 10.677655 6.201632 3.913725 2.661754 1.929642 1.474180
#> 12  83 37.71285 18.51298  9.908623 5.800306 3.698414 2.543635 1.863626 1.436691
#> 17  92 40.88788 19.62643 10.291830 5.926929 3.736897 2.553271 1.864739 1.435868
#> 31  76 35.26124 17.71995  9.696240 5.776095 3.724967 2.576328 1.890329 1.455856
#>         1.8         2
#> 21 1.176759 0.9686058
#> 35 0.932495 0.7983976
#> 33 1.158957 0.9586946
#> 5  1.175553 0.9678267
#> 16 1.175935 0.9686598
#> 7  1.173083 0.9672014
#> 14 1.182226 0.9718626
#> 22 1.153640 0.9548316
#> 48 1.174583 0.9676412
#> 12 1.152998 0.9550599
#> 17 1.152122 0.9545126
#> 31 1.165972 0.9635807
diversity(BCI[i,],"simpson") == x1[["2"]]
#>   21   35   33    5   16    7   14   22   48   12   17   31 
#> TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
x2 <- tsallis(BCI[i,],norm=TRUE)
x2
#>    0       0.2       0.4       0.6       0.8         1       1.2       1.4
#> 21 1 0.9070695 0.8546849 0.8354510 0.8415091 0.8639438 0.8936072 0.9229664
#> 35 1 0.7992527 0.6748192 0.6083778 0.5865659 0.5978626 0.6310844 0.6757387
#> 33 1 0.8944221 0.8348478 0.8119231 0.8167002 0.8397172 0.8714789 0.9040191
#> 5  1 0.9038026 0.8502581 0.8308757 0.8372445 0.8602030 0.8904867 0.9204832
#> 16 1 0.9070033 0.8547937 0.8356888 0.8417570 0.8641446 0.8938075 0.9232641
#> 7  1 0.9094899 0.8600339 0.8427936 0.8492653 0.8706729 0.8985570 0.9261475
#> 14 1 0.9134203 0.8653450 0.8483542 0.8548076 0.8762317 0.9040571 0.9313075
#> 22 1 0.8945746 0.8334045 0.8082215 0.8109045 0.8325271 0.8637842 0.8965990
#> 48 1 0.9106932 0.8602635 0.8414256 0.8466675 0.8676229 0.8957424 0.9239173
#> 12 1 0.8971858 0.8367628 0.8114291 0.8136186 0.8347024 0.8655269 0.8980648
#> 17 1 0.8945733 0.8308281 0.8026305 0.8032747 0.8244489 0.8566981 0.8913214
#> 31 1 0.9012670 0.8472634 0.8281644 0.8347227 0.8575355 0.8875792 0.9175864
#>          1.6       1.8         2
#> 21 0.9475331 0.9658643 0.9784896
#> 35 0.7233557 0.7684000 0.8081341
#> 33 0.9322682 0.9542083 0.9699733
#> 5  0.9456409 0.9644767 0.9775050
#> 16 0.9479928 0.9664763 0.9791887
#> 7  0.9493990 0.9669335 0.9791421
#> 14 0.9538903 0.9705576 0.9818818
#> 22 0.9256369 0.9486075 0.9654409
#> 48 0.9477911 0.9658285 0.9783927
#> 12 0.9269514 0.9498277 0.9665666
#> 17 0.9223038 0.9469046 0.9648877
#> 31 0.9431243 0.9625793 0.9762594
mod1 <- tsallisaccum(BCI[i,])
persp(mod1)

mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
```
