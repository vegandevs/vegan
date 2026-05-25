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
#> 30  96 42.46822 20.35888 10.670374 6.133357 3.851598 2.617135 1.900131 1.455342
#> 39  83 37.13183 17.98063  9.524369 5.542801 3.530494 2.435169 1.793695 1.391532
#> 49  90 40.07436 19.37264 10.256218 5.958648 3.778851 2.587992 1.889436 1.452206
#> 7   81 37.47840 18.73414 10.172688 6.004908 3.836811 2.631757 1.918091 1.469869
#> 29  85 38.51838 18.81994 10.009321 5.819419 3.688721 2.526664 1.847408 1.423613
#> 20  99 44.47873 21.56561 11.366723 6.530095 4.077327 2.746035 1.974163 1.498140
#> 41 101 45.01956 21.66615 11.353495 6.499317 4.052495 2.729799 1.964435 1.492566
#> 43  85 38.38228 18.77974 10.035449 5.867927 3.736254 2.565004 1.875593 1.443295
#> 40  79 34.49875 16.41508  8.624079 5.028357 3.234849 2.263379 1.692562 1.331193
#> 26  90 40.93055 20.10590 10.740144 6.250173 3.947749 2.684219 1.943942 1.483076
#> 34  91 40.76643 19.76149 10.454108 6.053163 3.821669 2.606302 1.896645 1.454634
#> 17  92 40.88788 19.62643 10.291830 5.926929 3.736897 2.553271 1.864739 1.435868
#>         1.8         2
#> 30 1.162753 0.9602659
#> 39 1.123737 0.9360204
#> 49 1.162505 0.9609552
#> 7  1.173083 0.9672014
#> 29 1.143249 0.9481041
#> 20 1.187661 0.9748589
#> 41 1.184546 0.9731442
#> 43 1.156578 0.9569632
#> 40 1.087270 0.9137131
#> 26 1.180039 0.9709567
#> 34 1.163015 0.9607876
#> 17 1.152122 0.9545126
diversity(BCI[i,],"simpson") == x1[["2"]]
#>   30   39   49    7   29   20   41   43   40   26   34   17 
#> TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
x2 <- tsallis(BCI[i,],norm=TRUE)
x2
#>    0       0.2       0.4       0.6       0.8         1       1.2       1.4
#> 30 1 0.8975551 0.8388546 0.8155960 0.8196218 0.8419326 0.8731640 0.9052902
#> 39 1 0.8833633 0.8127010 0.7799621 0.7774979 0.7968043 0.8286189 0.8643658
#> 49 1 0.8925852 0.8316222 0.8082153 0.8134945 0.8377231 0.8709199 0.9046665
#> 7  1 0.9094899 0.8600339 0.8427936 0.8492653 0.8706729 0.8985570 0.9261475
#> 29 1 0.8987625 0.8378269 0.8104423 0.8098004 0.8281171 0.8569302 0.8885453
#> 20 1 0.9168340 0.8714005 0.8563191 0.8638340 0.8853802 0.9124664 0.9383899
#> 41 1 0.9130269 0.8644349 0.8473152 0.8541311 0.8762204 0.9047063 0.9323803
#> 43 1 0.8955869 0.8360372 0.8125579 0.8165505 0.8387881 0.8699334 0.9021018
#> 40 1 0.8544019 0.7656946 0.7230720 0.7171847 0.7382083 0.7754969 0.8189349
#> 26 1 0.9116552 0.8630994 0.8463498 0.8532945 0.8751655 0.9033025 0.9307642
#> 34 1 0.8998773 0.8423760 0.8195118 0.8233665 0.8451677 0.8757762 0.9073385
#> 17 1 0.8945733 0.8308281 0.8026305 0.8032747 0.8244489 0.8566981 0.8913214
#>          1.6       1.8         2
#> 30 0.9331704 0.9547770 0.9702687
#> 39 0.8978144 0.9257229 0.9472977
#> 49 0.9336632 0.9558968 0.9716325
#> 7  0.9493990 0.9669335 0.9791421
#> 29 0.9175446 0.9412749 0.9592583
#> 20 0.9594192 0.9746101 0.9847059
#> 41 0.9550896 0.9716606 0.9827793
#> 43 0.9302296 0.9522490 0.9682216
#> 40 0.8608104 0.8967433 0.9252791
#> 26 0.9535106 0.9703141 0.9817451
#> 34 0.9347873 0.9560842 0.9713457
#> 17 0.9223038 0.9469046 0.9648877
mod1 <- tsallisaccum(BCI[i,])
persp(mod1)

mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
```
