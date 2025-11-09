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

Plotting methods and accumulation routines are based on functions
[`renyi`](https://vegandevs.github.io/vegan/reference/renyi.md) and
[`renyiaccum`](https://vegandevs.github.io/vegan/reference/renyi.md). An
object of class `tsallisaccum` can be displayed with dynamic 3D function
`rgl.renyiaccum` in the vegan3d package. See also settings for
[`persp`](https://rdrr.io/r/graphics/persp.html).

## Examples

``` r
data(BCI)
i <- sample(nrow(BCI), 12)
x1 <- tsallis(BCI[i,])
x1
#>      0      0.2      0.4       0.6      0.8        1      1.2      1.4      1.6
#> 27  98 43.84660 21.16041 11.110930 6.371476 3.980281 2.687154 1.938584 1.476662
#> 24  94 42.55654 20.75558 10.992840 6.343926 3.979427 2.692702 1.944478 1.481404
#> 40  79 34.49875 16.41508  8.624079 5.028357 3.234849 2.263379 1.692562 1.331193
#> 34  91 40.76643 19.76149 10.454108 6.053163 3.821669 2.606302 1.896645 1.454634
#> 26  90 40.93055 20.10590 10.740144 6.250173 3.947749 2.684219 1.943942 1.483076
#> 4   93 41.96845 20.49096 10.893385 6.315044 3.976563 2.697452 1.950144 1.486011
#> 20  99 44.47873 21.56561 11.366723 6.530095 4.077327 2.746035 1.974163 1.498140
#> 5  100 44.20598 21.17686 11.081257 6.350040 3.969940 2.683412 1.937940 1.477217
#> 32  87 39.46244 19.31897 10.289732 5.982060 3.784873 2.584344 1.882462 1.445176
#> 49  90 40.07436 19.37264 10.256218 5.958648 3.778851 2.587992 1.889436 1.452206
#> 17  92 40.88788 19.62643 10.291830 5.926929 3.736897 2.553271 1.864739 1.435868
#> 44  80 36.52805 18.09206  9.782209 5.776652 3.705016 2.555595 1.873754 1.443791
#>         1.8         2
#> 27 1.174682 0.9669962
#> 24 1.178152 0.9694268
#> 40 1.087270 0.9137131
#> 34 1.163015 0.9607876
#> 26 1.180039 0.9709567
#> 4  1.181429 0.9716117
#> 20 1.187661 0.9748589
#> 5  1.175553 0.9678267
#> 32 1.156659 0.9565267
#> 49 1.162505 0.9609552
#> 17 1.152122 0.9545126
#> 44 1.157572 0.9578733
diversity(BCI[i,],"simpson") == x1[["2"]]
#>   27   24   40   34   26    4   20    5   32   49   17   44 
#> TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
plot(x1)

x2 <- tsallis(BCI[i,],norm=TRUE)
x2
#>    0       0.2       0.4       0.6       0.8         1       1.2       1.4
#> 27 1 0.9112896 0.8605499 0.8410587 0.8456725 0.8661975 0.8940909 0.9221774
#> 24 1 0.9149445 0.8667059 0.8486326 0.8536824 0.8738548 0.9008891 0.9279024
#> 40 1 0.8544019 0.7656946 0.7230720 0.7171847 0.7382083 0.7754969 0.8189349
#> 34 1 0.8998773 0.8423760 0.8195118 0.8233665 0.8451677 0.8757762 0.9073385
#> 26 1 0.9116552 0.8630994 0.8463498 0.8532945 0.8751655 0.9033025 0.9307642
#> 4  1 0.9101799 0.8614882 0.8452139 0.8528120 0.8752597 0.9037665 0.9313684
#> 20 1 0.9168340 0.8714005 0.8563191 0.8638340 0.8853802 0.9124664 0.9383899
#> 5  1 0.9038026 0.8502581 0.8308757 0.8372445 0.8602030 0.8904867 0.9204832
#> 32 1 0.9035303 0.8473994 0.8239937 0.8259777 0.8453403 0.8737041 0.9037309
#> 49 1 0.8925852 0.8316222 0.8082153 0.8134945 0.8377231 0.8709199 0.9046665
#> 17 1 0.8945733 0.8308281 0.8026305 0.8032747 0.8244489 0.8566981 0.8913214
#> 44 1 0.8954397 0.8371683 0.8152612 0.8204162 0.8431126 0.8740718 0.9056624
#>          1.6       1.8         2
#> 27 0.9460502 0.9641597 0.9768635
#> 24 0.9507028 0.9678509 0.9797399
#> 40 0.8608104 0.8967433 0.9252791
#> 34 0.9347873 0.9560842 0.9713457
#> 26 0.9535106 0.9703141 0.9817451
#> 4  0.9540821 0.9707655 0.9820591
#> 20 0.9594192 0.9746101 0.9847059
#> 5  0.9456409 0.9644767 0.9775050
#> 32 0.9304962 0.9518103 0.9675213
#> 49 0.9336632 0.9558968 0.9716325
#> 17 0.9223038 0.9469046 0.9648877
#> 44 0.9330830 0.9544340 0.9698467
plot(x2)

mod1 <- tsallisaccum(BCI[i,])
plot(mod1, as.table=TRUE, col = c(1, 2, 2))

persp(mod1)

mod2 <- tsallisaccum(BCI[i,], norm=TRUE)
persp(mod2,theta=100,phi=30)
```
