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
#>  [1,]  3 162.3154 142.8448 184.6034 11.139508
#>  [2,]  4 175.0131 157.5132 198.7598 11.144021
#>  [3,]  5 183.0619 161.5945 201.4741 10.714205
#>  [4,]  6 189.0808 167.0091 210.9952 11.414785
#>  [5,]  7 193.3907 171.0579 214.4484 11.964836
#>  [6,]  8 197.9428 177.2255 223.3080 12.335888
#>  [7,]  9 201.5243 181.2740 227.1707 12.009374
#>  [8,] 10 203.5474 185.3626 223.9934 10.610820
#>  [9,] 11 207.4527 187.6362 233.7068 12.282477
#> [10,] 12 209.6998 190.1973 237.4029 11.865062
#> [11,] 13 212.3383 191.6459 238.3554 11.881198
#> [12,] 14 214.9288 193.3994 238.6777 12.423549
#> [13,] 15 216.7426 196.3209 242.2528 11.565688
#> [14,] 16 218.7221 196.6856 245.2822 11.927381
#> [15,] 17 220.6137 200.2000 248.3464 12.582336
#> [16,] 18 222.0628 200.8186 247.9958 12.699908
#> [17,] 19 223.0149 204.1625 248.0648 14.074630
#> [18,] 20 224.7878 207.7701 253.6955 13.465166
#> [19,] 21 226.0542 207.5765 249.3943 13.991082
#> [20,] 22 226.9677 209.6013 246.6202 11.288654
#> [21,] 23 227.3646 212.7728 248.8616 10.141657
#> [22,] 24 230.0540 213.5420 257.0580 11.528441
#> [23,] 25 231.6164 212.9304 257.8702 12.310503
#> [24,] 26 233.4008 214.0740 264.9190 14.194214
#> [25,] 27 235.0139 215.6205 261.6563 14.769726
#> [26,] 28 235.0408 215.2745 264.2844 14.154506
#> [27,] 29 235.0082 214.9416 264.7772 13.855116
#> [28,] 30 235.8563 217.3312 266.4591 13.166844
#> [29,] 31 236.1911 218.7744 263.9032 12.948101
#> [30,] 32 236.5633 219.0860 265.7853 12.235590
#> [31,] 33 236.7445 220.2854 263.9819 10.658291
#> [32,] 34 237.2423 220.9470 259.0892 10.517415
#> [33,] 35 237.3058 223.7120 258.0954  9.421569
#> [34,] 36 237.8443 223.8337 255.9057  9.440256
#> [35,] 37 236.9757 225.5187 255.4467  8.640692
#> [36,] 38 237.0225 225.5428 258.2641  8.773751
#> [37,] 39 237.5829 225.2564 260.7413  8.978335
#> [38,] 40 237.3777 225.7997 252.1612  7.779280
#> [39,] 41 237.5492 225.4600 255.2445  7.863511
#> [40,] 42 237.6753 224.5829 256.5710  8.452873
#> [41,] 43 237.2773 226.3343 254.3522  8.087743
#> [42,] 44 237.3751 226.0687 254.3531  7.701980
#> [43,] 45 237.2559 226.2559 255.9407  7.012555
#> [44,] 46 237.2378 227.2688 252.4475  6.346579
#> [45,] 47 237.3622 228.3434 252.4419  5.795966
#> [46,] 48 237.0146 228.7459 247.0639  4.338825
#> [47,] 49 236.4281 229.5569 245.4082  2.796530
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
