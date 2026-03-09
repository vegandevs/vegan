# Extrapolated Species Richness in a Species Pool

The functions estimate the extrapolated species richness in a species
pool, or the number of unobserved species. Function `specpool` is based
on incidences in sample sites, and gives a single estimate for a
collection of sample sites (matrix). Function `estimateR` is based on
abundances (counts) on single sample site.

## Usage

``` r
specpool(x, pool, smallsample = TRUE)
estimateR(x, ...)
specpool2vect(X, index = c("jack1","jack2", "chao", "boot","Species"))
estaccumR(x, permutations = 100, parallel = getOption("mc.cores"))
# S3 method for class 'poolaccum'
summary(object, display, alpha = 0.05, ...)
```

## Arguments

- x:

  Data frame or matrix with species data or the analysis result for
  `plot` function.

- pool:

  A vector giving a classification for pooling the sites in the species
  data. If missing, all sites are pooled together.

- smallsample:

  Use small sample correction \\(N-1)/N\\, where \\N\\ is the number of
  sites within the `pool`.

- X, object:

  A `specpool` result object.

- index:

  The selected index of extrapolated richness.

- permutations:

  Usually an integer giving the number permutations, but can also be a
  list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or a
  permutation matrix where each row gives the permuted indices.

- parallel:

  Number of parallel processes or a predefined socket cluster. With
  `parallel = 1` uses ordinary, non-parallel processing. The parallel
  processing is done with parallel package.

- display:

  Indices to be displayed.

- alpha:

  Level of quantiles shown. This proportion will be left outside
  symmetric limits.

- ...:

  Other parameters (not used).

## Details

Many species will always remain unseen or undetected in a collection of
sample plots. The function uses some popular ways of estimating the
number of these unseen species and adding them to the observed species
richness (Palmer 1990, Colwell & Coddington 1994).

The incidence-based estimates in `specpool` use the frequencies of
species in a collection of sites. In the following, \\S_P\\ is the
extrapolated richness in a pool, \\S_0\\ is the observed number of
species in the collection, \\a_1\\ and \\a_2\\ are the number of species
occurring only in one or only in two sites in the collection, \\p_i\\ is
the frequency of species \\i\\, and \\N\\ is the number of sites in the
collection. The variants of extrapolated richness in `specpool` are:

|                        |                                                                    |
|------------------------|--------------------------------------------------------------------|
| Chao                   | \\S_P = S_0 + \frac{a_1^2}{2 a_2}\frac{N-1}{N}\\                   |
| Chao bias-corrected    | \\S_P = S_0 + \frac{a_1(a_1-1)}{2(a_2+1)} \frac{N-1}{N}\\          |
| First order jackknife  | \\S_P = S_0 + a_1 \frac{N-1}{N}\\                                  |
| Second order jackknife | \\S_P = S_0 + a_1 \frac{2N - 3}{N} - a_2 \frac{(N-2)^2}{N (N-1)}\\ |
| Bootstrap              | \\S_P = S_0 + \sum\_{i=1}^{S_0} (1 - p_i)^N\\                      |

`specpool` normally uses basic Chao equation, but when there are no
doubletons (\\a2=0\\) it switches to bias-corrected version. In that
case the Chao equation simplifies to \\S_0 + \frac{1}{2} a_1 (a_1-1)
\frac{N-1}{N}\\.

The abundance-based estimates in `estimateR` use counts (numbers of
individuals) of species in a single site. If called for a matrix or data
frame, the function will give separate estimates for each site. The two
variants of extrapolated richness in `estimateR` are bias-corrected Chao
and ACE (O'Hara 2005, Chiu et al. 2014). The Chao estimate is similar as
the bias corrected one above, but \\a_i\\ refers to the number of
species with abundance \\i\\ instead of number of sites, and the
small-sample correction is not used. The ACE estimate is defined as:

|       |                                                                                                                                  |
|-------|----------------------------------------------------------------------------------------------------------------------------------|
| ACE   | \\S_P = S\_{abund} + \frac{S\_{rare}}{C\_{ace}}+ \frac{a_1}{C\_{ace}} \gamma^2\_{ace}\\                                          |
| where | \\C\_{ace} = 1 - \frac{a_1}{N\_{rare}}\\                                                                                         |
|       | \\\gamma^2\_{ace} = \max \left\[ \frac{S\_{rare} \sum\_{i=1}^{10} i(i-1)a_i}{C\_{ace} N\_{rare} (N\_{rare} - 1)}-1, 0 \right\]\\ |

Here \\a_i\\ refers to number of species with abundance \\i\\ and
\\S\_{rare}\\ is the number of rare species, \\S\_{abund}\\ is the
number of abundant species, with an arbitrary threshold of abundance 10
for rare species, and \\N\_{rare}\\ is the number of individuals in rare
species.

Functions estimate the standard errors of the estimates. These only
concern the number of added species, and assume that there is no
variance in the observed richness. The equations of standard errors are
too complicated to be reproduced in this help page, but they can be
studied in the R source code of the function and are discussed in the
[`vignette`](https://rdrr.io/r/utils/vignette.html) that can be read
with the `browseVignettes("vegan")`. The standard error are based on the
following sources: Chiu et al. (2014) for the Chao estimates and Smith
and van Belle (1984) for the first-order Jackknife and the bootstrap
(second-order jackknife is still missing). For the variance estimator of
\\S\_{ace}\\ see O'Hara (2005).

Functions `poolaccum` and `estaccumR` are similar to
[`specaccum`](https://vegandevs.github.io/vegan/reference/specaccum.md),
but estimate extrapolated richness indices of `specpool` or `estimateR`
in addition to number of species for random ordering of sampling units.
Function `specpool` uses presence data and `estaccumR` count data. The
functions share `summary` and `plot` methods. The `summary` returns
quantile envelopes of permutations corresponding the given level of
`alpha` and standard deviation of permutations for each sample size.
NB., these are not based on standard deviations estimated within
`specpool` or `estimateR`, but they are based on permutations.

[ggvegan](https://CRAN.R-project.org/package=ggvegan) provides function
`autoplot` to show the results of `poolaccum`.

## Value

Function `specpool` returns a data frame with entries for observed
richness and each of the indices for each class in `pool` vector. The
utility function `specpool2vect` maps the pooled values into a vector
giving the value of selected `index` for each original site. Function
`estimateR` returns the estimates and their standard errors for each
site. Functions `poolaccum` and `estimateR` return matrices of
permutation results for each richness estimator, the vector of sample
sizes and a table of `means` of permutations for each estimator.

## References

Chao, A. (1987). Estimating the population size for capture-recapture
data with unequal catchability. *Biometrics* 43, 783–791.

Chiu, C.H., Wang, Y.T., Walther, B.A. & Chao, A. (2014). Improved
nonparametric lower bound of species richness via a modified Good-Turing
frequency formula. *Biometrics* 70, 671–682.

Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial
biodiversity through extrapolation. *Phil. Trans. Roy. Soc. London* B
345, 101–118.

O'Hara, R.B. (2005). Species richness estimators: how many species can
dance on the head of a pin? *J. Anim. Ecol.* 74, 375–386.

Palmer, M.W. (1990). The estimation of species richness by
extrapolation. *Ecology* 71, 1195–1198.

Smith, E.P & van Belle, G. (1984). Nonparametric estimation of species
richness. *Biometrics* 40, 119–129.

## Author

Bob O'Hara (`estimateR`) and Jari Oksanen.

## Note

The functions are based on assumption that there is a species pool: The
community is closed so that there is a fixed pool size \\S_P\\. In
general, the functions give only the lower limit of species richness:
the real richness is \\S \>= S_P\\, and there is a consistent bias in
the estimates. Even the bias-correction in Chao only reduces the bias,
but does not remove it completely (Chiu et al. 2014).

Optional small sample correction was added to `specpool` in vegan 2.2-0.
It was not used in the older literature (Chao 1987), but it is
recommended recently (Chiu et al. 2014).

## See also

[`veiledspec`](https://vegandevs.github.io/vegan/reference/fisherfit.md),
[`diversity`](https://vegandevs.github.io/vegan/reference/diversity.md),
[`beals`](https://vegandevs.github.io/vegan/reference/beals.md),
[`specaccum`](https://vegandevs.github.io/vegan/reference/specaccum.md).

## Examples

``` r
data(dune)
data(dune.env)
pool <- with(dune.env, specpool(dune, Management))
pool
#>    Species     chao   chao.se    jack1 jack1.se    jack2     boot  boot.se n
#> BF      16 17.19048 1.5895675 19.33333 2.211083 19.83333 17.74074 1.646379 3
#> HF      21 21.51429 0.9511693 23.40000 1.876166 22.05000 22.56864 1.821518 5
#> NM      21 22.87500 2.1582871 26.00000 3.291403 25.73333 23.77696 2.300982 6
#> SF      21 29.88889 8.6447967 27.66667 3.496029 31.40000 23.99496 1.850288 6
op <- par(mfrow=c(1,2))
boxplot(specnumber(dune) ~ Management, data = dune.env,
        col = "hotpink", border = "cyan3")
boxplot(specnumber(dune)/specpool2vect(pool) ~ Management,
        data = dune.env, col = "hotpink", border = "cyan3")

par(op)
data(BCI)
## Accumulation model
pool <- poolaccum(BCI)
summary(pool, display = "chao")
#> $chao
#>        N     Chao     2.5%    97.5%   Std.Dev
#>  [1,]  3 161.9126 137.9651 181.8037 11.392972
#>  [2,]  4 174.5347 157.5132 197.4215 11.623127
#>  [3,]  5 181.6076 160.6902 202.5886 11.835495
#>  [4,]  6 187.8478 167.0091 209.0661 11.496356
#>  [5,]  7 192.2664 171.0579 215.0067 12.172911
#>  [6,]  8 196.6568 177.2255 220.4803 11.795564
#>  [7,]  9 201.1551 182.1626 227.5839 12.903996
#>  [8,] 10 203.6162 185.3626 226.1050 11.288463
#>  [9,] 11 207.5008 187.6362 233.3426 12.420970
#> [10,] 12 209.3699 190.1973 234.0910 11.556014
#> [11,] 13 212.2943 192.9427 239.2166 11.646958
#> [12,] 14 214.1672 194.4139 238.6777 11.851882
#> [13,] 15 216.4289 196.3209 244.4969 11.662370
#> [14,] 16 218.1599 196.6856 243.6845 12.223573
#> [15,] 17 220.0908 200.2000 249.0297 13.151683
#> [16,] 18 221.5332 200.8186 254.9685 13.846710
#> [17,] 19 222.4574 203.0215 249.0148 14.776242
#> [18,] 20 224.3713 205.8278 256.2941 14.017309
#> [19,] 21 225.5285 206.0792 256.9636 14.378280
#> [20,] 22 226.1096 208.0329 252.2394 11.948046
#> [21,] 23 226.6964 210.8317 250.9183 11.248721
#> [22,] 24 228.9428 212.4638 258.4252 11.848765
#> [23,] 25 230.3934 212.7523 260.7472 12.230439
#> [24,] 26 232.3091 213.9838 264.9190 14.121262
#> [25,] 27 233.5665 213.8450 261.6563 14.255698
#> [26,] 28 233.3065 214.7239 263.7886 13.267956
#> [27,] 29 233.4431 214.2811 262.7931 13.265664
#> [28,] 30 234.3525 214.8768 263.7592 13.091248
#> [29,] 31 234.9381 216.4895 263.2367 12.340616
#> [30,] 32 235.6385 217.8520 259.6366 11.414028
#> [31,] 33 236.3385 218.2835 255.1406 10.238112
#> [32,] 34 237.0426 219.3071 257.9771 10.066652
#> [33,] 35 237.5603 221.8241 254.4917  9.124579
#> [34,] 36 238.1143 222.1373 255.9057  9.391886
#> [35,] 37 237.5235 224.0933 255.4467  8.618626
#> [36,] 38 237.5961 225.7349 256.6979  8.489668
#> [37,] 39 237.9947 225.3127 259.5813  8.485005
#> [38,] 40 237.8337 225.7997 252.1065  7.602304
#> [39,] 41 237.6682 225.4600 253.0163  7.814481
#> [40,] 42 237.7471 224.5829 253.9004  8.489582
#> [41,] 43 237.2841 225.6792 254.6499  8.268247
#> [42,] 44 237.5236 225.8592 254.3531  7.855885
#> [43,] 45 237.4549 226.4246 256.6274  7.368402
#> [44,] 46 237.6986 227.2688 254.4983  6.424953
#> [45,] 47 237.3845 228.1245 252.4419  5.657572
#> [46,] 48 237.0711 228.7459 247.0639  4.475686
#> [47,] 49 236.4252 229.5569 242.9539  2.757967
#> [48,] 50 236.3732 236.3732 236.3732  0.000000
#> 
#> attr(,"class")
#> [1] "summary.poolaccum"
## Quantitative model
estimateR(BCI[1:5,])
#>                   1          2          3          4          5
#> S.obs     93.000000  84.000000  90.000000  94.000000 101.000000
#> S.chao1  117.473684 117.214286 141.230769 111.550000 136.000000
#> se.chao1  11.583785  15.918953  23.001405   8.919663  15.467344
#> S.ACE    122.848959 117.317307 134.669844 118.729941 137.114088
#> se.ACE     5.736054   5.571998   6.191618   5.367571   5.848474
```
