# Species Accumulation Curves

Function `specaccum` finds species accumulation curves or the number of
species for a certain number of sampled sites or individuals.

## Usage

``` r
specaccum(comm, method = "exact", permutations = 100,
          conditioned =TRUE, gamma = "jack1",  w = NULL, subset, ...)
# S3 method for class 'specaccum'
plot(x, add = FALSE, random = FALSE, ci = 2, 
    ci.type = c("bar", "line", "polygon"), col = par("fg"), lty = 1,
    ci.col = col, ci.lty = 1, ci.length = 0, xlab, ylab = x$method, ylim,
    xvar = c("sites", "individuals", "effort"), ...)
# S3 method for class 'specaccum'
boxplot(x, add = FALSE, ...)
fitspecaccum(object, model, method = "random", ...)
# S3 method for class 'fitspecaccum'
plot(x, col = par("fg"), lty = 1, xlab = "Sites", 
    ylab = x$method, ...) 
# S3 method for class 'specaccum'
predict(object, newdata, interpolation = c("linear", "spline"), ...)
# S3 method for class 'fitspecaccum'
predict(object, newdata, ...)
specslope(object, at)
```

## Arguments

- comm:

  Community data set.

- method:

  Species accumulation method (partial match). Method `"collector"` adds
  sites in the order they happen to be in the data, `"random"` adds
  sites in random order, `"exact"` finds the expected (mean) species
  richness, `"coleman"` finds the expected richness following Coleman et
  al. 1982, and `"rarefaction"` finds the mean when accumulating
  individuals instead of sites.

- permutations:

  Number of permutations with `method = "random"`. Usually an integer
  giving the number permutations, but can also be a list of control
  values for the permutations as returned by the function
  [`how`](https://rdrr.io/pkg/permute/man/how.html), or a permutation
  matrix where each row gives the permuted indices.

- conditioned:

  Estimation of standard deviation is conditional on the empirical
  dataset for the exact SAC

- gamma:

  Method for estimating the total extrapolated number of species in the
  survey area by function
  [`specpool`](https://vegandevs.github.io/vegan/reference/specpool.md)

- w:

  Weights giving the sampling effort.

- subset:

  logical expression indicating sites (rows) to keep: missing values are
  taken as `FALSE`.

- x:

  A `specaccum` result object

- add:

  Add to an existing graph.

- random:

  Draw each random simulation separately instead of drawing their
  average and confidence intervals.

- ci:

  Multiplier used to get confidence intervals from standard deviation
  (standard error of the estimate). Value `ci = 0` suppresses drawing
  confidence intervals.

- ci.type:

  Type of confidence intervals in the graph: `"bar"` draws vertical
  bars, `"line"` draws lines, and `"polygon"` draws a shaded area.

- col:

  Colour for drawing lines.

- lty:

  line type (see [`par`](https://rdrr.io/r/graphics/par.html)).

- ci.col:

  Colour for drawing lines or filling the `"polygon"`.

- ci.lty:

  Line type for confidence intervals or border of the `"polygon"`.

- ci.length:

  Length of horizontal bars (in inches) at the end of vertical bars with
  `ci.type = "bar"`.

- xlab,ylab:

  Labels for `x` (defaults `xvar`) and `y` axis.

- ylim:

  the y limits of the plot.

- xvar:

  Variable used for the horizontal axis: `"individuals"` can be used
  only with `method = "rarefaction"`.

- object:

  Either a community data set or fitted `specaccum` model.

- model:

  Nonlinear regression model
  ([`nls`](https://rdrr.io/r/stats/nls.html)). See Details.

- newdata:

  Optional data used in prediction interpreted as number of sampling
  units (sites). If missing, fitted values are returned.

- interpolation:

  Interpolation method used with `newdata`.

- at:

  Number of plots where the slope is evaluated. Can be a real number.

- ...:

  Other parameters to functions.

## Details

Species accumulation curves (SAC) are used to compare diversity
properties of community data sets using different accumulator functions.
The classic method is `"random"` which finds the mean SAC and its
standard deviation from random permutations of the data, or subsampling
without replacement (Gotelli & Colwell 2001). The `"exact"` method finds
the expected SAC using sample-based rarefaction method that has been
independently developed numerous times (Chiarucci et al. 2008) and it is
often known as Mao Tau estimate (Colwell et al. 2012). The unconditional
standard deviation for the exact SAC represents a moment-based
estimation that is not conditioned on the empirical data set (sd for all
samples \> 0). The unconditional standard deviation is based on an
estimation of the extrapolated number of species in the survey area
(a.k.a. gamma diversity), as estimated by function
[`specpool`](https://vegandevs.github.io/vegan/reference/specpool.md).
The conditional standard deviation that was developed by Jari Oksanen
(not published, sd=0 for all samples). Method `"coleman"` finds the
expected SAC and its standard deviation following Coleman et al. (1982).
All these methods are based on sampling sites without replacement. In
contrast, the `method = "rarefaction"` finds the expected species
richness and its standard deviation by sampling individuals instead of
sites. It achieves this by applying function
[`rarefy`](https://vegandevs.github.io/vegan/reference/rarefy.md) with
number of individuals corresponding to average number of individuals per
site.

Methods `"random"` and `"collector"` can take weights (`w`) that give
the sampling effort for each site. The weights `w` do not influence the
order the sites are accumulated, but only the value of the sampling
effort so that not all sites are equal. The summary results are
expressed against sites even when the accumulation uses weights (methods
`"random"`, `"collector"`), or is based on individuals
(`"rarefaction"`). The actual sampling effort is given as item `Effort`
or `Individuals` in the printed result. For weighted `"random"` method
the effort refers to the average effort per site, or sum of weights per
number of sites. With weighted `method = "random"`, the averaged species
richness is found from linear interpolation of single random
permutations. Therefore at least the first value (and often several
first) have `NA` richness, because these values cannot be interpolated
in all cases but should be extrapolated. The `plot` function defaults to
display the results as scaled to sites, but this can be changed
selecting `xvar = "effort"` (weighted methods) or `xvar = "individuals"`
(with `method = "rarefaction"`).

The `summary` and `boxplot` methods are available for
`method = "random"`.

Function `predict` for `specaccum` can return the values corresponding
to `newdata`. With `method` `"exact"`, `"rarefaction"` and `"coleman"`
the function uses analytic equations for interpolated non-integer
values, and for other methods linear
([`approx`](https://rdrr.io/r/stats/approxfun.html)) or spline
([`spline`](https://rdrr.io/r/stats/splinefun.html)) interpolation. If
`newdata` is not given, the function returns the values corresponding to
the data. NB., the fitted values with `method="rarefaction"` are based
on rounded integer counts, but `predict` can use fractional non-integer
counts with `newdata` and give slightly different results.

Function `fitspecaccum` fits a nonlinear
([`nls`](https://rdrr.io/r/stats/nls.html)) self-starting species
accumulation model. The input `object` can be a result of `specaccum` or
a community in data frame. In the latter case the function first fits a
`specaccum` model and then proceeds with fitting the nonlinear model.
The function can apply a limited set of nonlinear regression models
suggested for species-area relationship (Dengler 2009). All these are
[`selfStart`](https://rdrr.io/r/stats/selfStart.html) models. The
permissible alternatives are `"arrhenius"`
([`SSarrhenius`](https://vegandevs.github.io/vegan/reference/SSarrhenius.md)),
`"gleason"`
([`SSgleason`](https://vegandevs.github.io/vegan/reference/SSarrhenius.md)),
`"gitay"`
([`SSgitay`](https://vegandevs.github.io/vegan/reference/SSarrhenius.md)),
`"lomolino"`
([`SSlomolino`](https://vegandevs.github.io/vegan/reference/SSarrhenius.md))
of vegan package. In addition the following standard R models are
available: `"asymp"`
([`SSasymp`](https://rdrr.io/r/stats/SSasymp.html)), `"gompertz"`
([`SSgompertz`](https://rdrr.io/r/stats/SSgompertz.html)),
`"michaelis-menten"`
([`SSmicmen`](https://rdrr.io/r/stats/SSmicmen.html)), `"logis"`
([`SSlogis`](https://rdrr.io/r/stats/SSlogis.html)), `"weibull"`
([`SSweibull`](https://rdrr.io/r/stats/SSweibull.html)). See these
functions for model specification and details.

When weights `w` were used the fit is based on accumulated effort and in
`model = "rarefaction"` on accumulated number of individuals. The `plot`
is still based on sites, unless other alternative is selected with
`xvar`.

Function `predict` for `fitspecaccum` uses
[`predict.nls`](https://rdrr.io/r/stats/predict.nls.html), and you can
pass all arguments to that function. In addition, `fitted`, `residuals`,
`nobs`, `coef`, `AIC`, `logLik` and `deviance` work on the result
object.

Function `specslope` evaluates the derivative of the species
accumulation curve at given number of sample plots, and gives the rate
of increase in the number of species. The function works with
`specaccum` result object when this is based on analytic models
`"exact"`, `"rarefaction"` or `"coleman"`, and with non-linear
regression results of `fitspecaccum`.

Nonlinear regression may fail for any reason, and some of the
`fitspecaccum` models are fragile and may not succeed.

## Value

Function `specaccum` returns an object of class `"specaccum"`, and
`fitspecaccum` a model of class `"fitspecaccum"` that adds a few items
to the `"specaccum"` (see the end of the list below):

- call :

  Function call.

- method:

  Accumulator method.

- sites:

  Number of sites. For `method = "rarefaction"` this is the number of
  sites corresponding to a certain number of individuals and generally
  not an integer, and the average number of individuals is also returned
  in item `individuals`.

- effort:

  Average sum of weights corresponding to the number of sites when model
  was fitted with argument `w`

- richness:

  The number of species corresponding to number of sites. With
  `method = "collector"` this is the observed richness, for other
  methods the average or expected richness.

- sd:

  The standard deviation of SAC (or its standard error). This is `NULL`
  in `method = "collector"`, and it is estimated from permutations in
  `method = "random"`, and from analytic equations in other methods.

- perm:

  Permutation results with `method = "random"` and `NULL` in other
  cases. Each column in `perm` holds one permutation.

- weights:

  Matrix of accumulated weights corresponding to the columns of the
  `perm` matrix when model was fitted with argument `w`.

- fitted, residuals, coefficients:

  Only in `fitspecacum`: fitted values, residuals and nonlinear model
  coefficients. For `method = "random"` these are matrices with a column
  for each random accumulation.

- models:

  Only in `fitspecaccum`: list of fitted
  [`nls`](https://rdrr.io/r/stats/nls.html) models (see Examples on
  accessing these models).

## References

Chiarucci, A., Bacaro, G., Rocchini, D. & Fattorini, L. (2008).
Discovering and rediscovering the sample-based rarefaction formula in
the ecological literature. *Commun. Ecol.* 9: 121–123.

Coleman, B.D, Mares, M.A., Willis, M.R. & Hsieh, Y. (1982). Randomness,
area and species richness. *Ecology* 63: 1121–1133.

Colwell, R.K., Chao, A., Gotelli, N.J., Lin, S.Y., Mao, C.X., Chazdon,
R.L. & Longino, J.T. (2012). Models and estimators linking
individual-based and sample-based rarefaction, extrapolation and
comparison of assemblages. *J. Plant Ecol.* 5: 3–21.

Dengler, J. (2009). Which function describes the species-area
relationship best? A review and empirical evaluation. *Journal of
Biogeography* 36, 728–744.

Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
procedures and pitfalls in measurement and comparison of species
richness. *Ecol. Lett.* 4, 379–391.

## Author

Roeland Kindt <r.kindt@cgiar.org> and Jari Oksanen.

## Note

The SAC with `method = "exact"` was developed by Roeland Kindt, and its
standard deviation by Jari Oksanen (both are unpublished). The
`method = "coleman"` underestimates the SAC because it does not handle
properly sampling without replacement. Further, its standard deviation
does not take into account species correlations, and is generally too
low.

## See also

[`rarefy`](https://vegandevs.github.io/vegan/reference/rarefy.md) and
[`rrarefy`](https://vegandevs.github.io/vegan/reference/rarefy.md) are
related individual based models. Other accumulation models are
[`poolaccum`](https://vegandevs.github.io/vegan/reference/specpool.md)
for extrapolated richness, and
[`renyiaccum`](https://vegandevs.github.io/vegan/reference/renyi.md) and
[`tsallisaccum`](https://vegandevs.github.io/vegan/reference/tsallis.md)
for diversity indices. Underlying graphical functions are
[`boxplot`](https://rdrr.io/r/graphics/boxplot.html),
[`matlines`](https://rdrr.io/r/graphics/matplot.html),
[`segments`](https://rdrr.io/r/graphics/segments.html) and
[`polygon`](https://rdrr.io/r/graphics/polygon.html).

## Examples

``` r
data(BCI)
sp1 <- specaccum(BCI)
#> Warning: the standard deviation is zero
sp2 <- specaccum(BCI, "random")
sp2
#> Species Accumulation Curve
#> Accumulation method: random, with 100 permutations
#> Call: specaccum(comm = BCI, method = "random") 
#> 
#>                                                                              
#> Sites     1.00000   2.00000   3.00000   4.00000   5.00000   6.00000   7.00000
#> Richness 90.03000 121.37000 139.25000 150.84000 158.97000 165.33000 170.98000
#> sd        6.24978   6.45521   5.58384   5.44897   5.42805   5.18556   5.08907
#>                                                                               
#> Sites      8.00000   9.00000  10.00000  11.00000  12.00000  13.00000  14.00000
#> Richness 175.73000 179.67000 182.51000 185.38000 188.08000 190.71000 192.84000
#> sd         5.05696   4.73767   4.44835   4.41938   4.22518   3.99316   3.88657
#>                                                                            
#> Sites     15.00000  16.00000  17.0000  18.00000  19.000  20.00000  21.00000
#> Richness 194.61000 196.34000 197.8000 199.30000 200.560 202.04000 203.29000
#> sd         4.11917   4.25505   4.1827   4.05642   4.036   3.92073   3.78005
#>                                                                              
#> Sites     22.00000  23.00000  24.00000  25.0000  26.00000  27.00000  28.00000
#> Richness 204.47000 205.50000 206.75000 207.7600 208.79000 209.95000 210.74000
#> sd         3.73748   3.69411   3.53446   3.4848   3.46788   3.42709   3.50647
#>                                                                               
#> Sites     29.00000  30.00000  31.00000  32.00000  33.00000  34.00000  35.00000
#> Richness 211.72000 212.54000 213.49000 214.29000 215.28000 215.93000 216.76000
#> sd         3.43505   3.37076   3.36799   3.09543   3.13688   3.07894   2.73444
#>                                                                               
#> Sites     36.00000  37.00000  38.00000  39.00000  40.00000  41.00000  42.00000
#> Richness 217.36000 217.95000 218.56000 219.21000 219.98000 220.55000 220.97000
#> sd         2.61897   2.44278   2.36267   2.19869   2.08399   1.98161   2.01737
#>                                                                             
#> Sites     43.0000  44.0000  45.00000  46.00000  47.00000  48.00000  49.00000
#> Richness 221.5100 222.0400 222.51000 223.08000 223.62000 224.12000 224.53000
#> sd         1.9201   1.9065   1.78939   1.53531   1.30097   1.10353   0.82211
#>             
#> Sites     50
#> Richness 225
#> sd         0
summary(sp2)
#>  1 sites          2 sites         3 sites         4 sites         5 sites      
#>  Min.   : 77.00   Min.   :110.0   Min.   :126.0   Min.   :140.0   Min.   :144  
#>  1st Qu.: 84.75   1st Qu.:116.0   1st Qu.:136.0   1st Qu.:146.0   1st Qu.:156  
#>  Median : 90.00   Median :121.5   Median :139.5   Median :151.0   Median :158  
#>  Mean   : 90.03   Mean   :121.4   Mean   :139.2   Mean   :150.8   Mean   :159  
#>  3rd Qu.: 93.00   3rd Qu.:126.2   3rd Qu.:143.0   3rd Qu.:155.0   3rd Qu.:163  
#>  Max.   :105.00   Max.   :135.0   Max.   :152.0   Max.   :162.0   Max.   :172  
#>  6 sites         7 sites         8 sites         9 sites        
#>  Min.   :151.0   Min.   :158.0   Min.   :161.0   Min.   :168.0  
#>  1st Qu.:162.0   1st Qu.:168.0   1st Qu.:172.8   1st Qu.:177.0  
#>  Median :165.0   Median :171.0   Median :175.0   Median :179.0  
#>  Mean   :165.3   Mean   :171.0   Mean   :175.7   Mean   :179.7  
#>  3rd Qu.:169.0   3rd Qu.:174.2   3rd Qu.:179.0   3rd Qu.:182.0  
#>  Max.   :179.0   Max.   :182.0   Max.   :188.0   Max.   :192.0  
#>  10 sites        11 sites        12 sites        13 sites       
#>  Min.   :170.0   Min.   :176.0   Min.   :179.0   Min.   :181.0  
#>  1st Qu.:180.0   1st Qu.:182.8   1st Qu.:185.0   1st Qu.:188.0  
#>  Median :182.0   Median :185.0   Median :188.0   Median :191.0  
#>  Mean   :182.5   Mean   :185.4   Mean   :188.1   Mean   :190.7  
#>  3rd Qu.:185.0   3rd Qu.:188.0   3rd Qu.:191.2   3rd Qu.:193.0  
#>  Max.   :194.0   Max.   :195.0   Max.   :197.0   Max.   :200.0  
#>  14 sites        15 sites        16 sites        17 sites       
#>  Min.   :183.0   Min.   :184.0   Min.   :185.0   Min.   :186.0  
#>  1st Qu.:190.8   1st Qu.:192.0   1st Qu.:194.0   1st Qu.:196.0  
#>  Median :193.0   Median :194.5   Median :196.0   Median :198.0  
#>  Mean   :192.8   Mean   :194.6   Mean   :196.3   Mean   :197.8  
#>  3rd Qu.:195.0   3rd Qu.:197.0   3rd Qu.:198.2   3rd Qu.:200.0  
#>  Max.   :201.0   Max.   :204.0   Max.   :208.0   Max.   :210.0  
#>  18 sites        19 sites        20 sites      21 sites        22 sites       
#>  Min.   :189.0   Min.   :190.0   Min.   :191   Min.   :192.0   Min.   :194.0  
#>  1st Qu.:197.0   1st Qu.:198.0   1st Qu.:199   1st Qu.:201.0   1st Qu.:202.0  
#>  Median :199.0   Median :200.0   Median :202   Median :203.0   Median :204.0  
#>  Mean   :199.3   Mean   :200.6   Mean   :202   Mean   :203.3   Mean   :204.5  
#>  3rd Qu.:202.0   3rd Qu.:202.2   3rd Qu.:205   3rd Qu.:206.0   3rd Qu.:207.0  
#>  Max.   :210.0   Max.   :210.0   Max.   :210   Max.   :212.0   Max.   :212.0  
#>  23 sites        24 sites        25 sites        26 sites       
#>  Min.   :196.0   Min.   :197.0   Min.   :199.0   Min.   :199.0  
#>  1st Qu.:203.0   1st Qu.:204.8   1st Qu.:205.0   1st Qu.:207.0  
#>  Median :205.0   Median :207.0   Median :208.0   Median :209.0  
#>  Mean   :205.5   Mean   :206.8   Mean   :207.8   Mean   :208.8  
#>  3rd Qu.:208.0   3rd Qu.:209.0   3rd Qu.:210.0   3rd Qu.:211.0  
#>  Max.   :214.0   Max.   :215.0   Max.   :216.0   Max.   :216.0  
#>  27 sites        28 sites        29 sites        30 sites       
#>  Min.   :201.0   Min.   :203.0   Min.   :203.0   Min.   :203.0  
#>  1st Qu.:208.0   1st Qu.:208.0   1st Qu.:209.0   1st Qu.:210.0  
#>  Median :210.0   Median :211.0   Median :213.0   Median :213.0  
#>  Mean   :209.9   Mean   :210.7   Mean   :211.7   Mean   :212.5  
#>  3rd Qu.:212.0   3rd Qu.:213.0   3rd Qu.:214.0   3rd Qu.:215.0  
#>  Max.   :217.0   Max.   :217.0   Max.   :218.0   Max.   :218.0  
#>  31 sites        32 sites        33 sites        34 sites       
#>  Min.   :204.0   Min.   :208.0   Min.   :208.0   Min.   :209.0  
#>  1st Qu.:212.0   1st Qu.:212.0   1st Qu.:213.0   1st Qu.:214.0  
#>  Median :214.0   Median :214.5   Median :216.0   Median :216.0  
#>  Mean   :213.5   Mean   :214.3   Mean   :215.3   Mean   :215.9  
#>  3rd Qu.:216.0   3rd Qu.:216.0   3rd Qu.:217.2   3rd Qu.:218.0  
#>  Max.   :221.0   Max.   :221.0   Max.   :221.0   Max.   :222.0  
#>  35 sites        36 sites        37 sites        38 sites       
#>  Min.   :211.0   Min.   :211.0   Min.   :213.0   Min.   :214.0  
#>  1st Qu.:215.0   1st Qu.:216.0   1st Qu.:216.0   1st Qu.:217.0  
#>  Median :217.0   Median :217.0   Median :218.0   Median :218.5  
#>  Mean   :216.8   Mean   :217.4   Mean   :217.9   Mean   :218.6  
#>  3rd Qu.:219.0   3rd Qu.:219.0   3rd Qu.:219.0   3rd Qu.:220.0  
#>  Max.   :222.0   Max.   :222.0   Max.   :224.0   Max.   :224.0  
#>  39 sites        40 sites        41 sites        42 sites      43 sites       
#>  Min.   :215.0   Min.   :215.0   Min.   :216.0   Min.   :216   Min.   :216.0  
#>  1st Qu.:217.0   1st Qu.:218.8   1st Qu.:219.0   1st Qu.:220   1st Qu.:220.0  
#>  Median :219.0   Median :220.0   Median :221.0   Median :221   Median :222.0  
#>  Mean   :219.2   Mean   :220.0   Mean   :220.6   Mean   :221   Mean   :221.5  
#>  3rd Qu.:221.0   3rd Qu.:222.0   3rd Qu.:222.0   3rd Qu.:222   3rd Qu.:223.0  
#>  Max.   :224.0   Max.   :225.0   Max.   :225.0   Max.   :225   Max.   :225.0  
#>  44 sites      45 sites        46 sites        47 sites        48 sites       
#>  Min.   :216   Min.   :217.0   Min.   :218.0   Min.   :220.0   Min.   :220.0  
#>  1st Qu.:221   1st Qu.:221.0   1st Qu.:222.0   1st Qu.:223.0   1st Qu.:224.0  
#>  Median :222   Median :223.0   Median :223.0   Median :224.0   Median :224.0  
#>  Mean   :222   Mean   :222.5   Mean   :223.1   Mean   :223.6   Mean   :224.1  
#>  3rd Qu.:223   3rd Qu.:224.0   3rd Qu.:224.0   3rd Qu.:225.0   3rd Qu.:225.0  
#>  Max.   :225   Max.   :225.0   Max.   :225.0   Max.   :225.0   Max.   :225.0  
#>  49 sites        50 sites     
#>  Min.   :222.0   Min.   :225  
#>  1st Qu.:224.0   1st Qu.:225  
#>  Median :225.0   Median :225  
#>  Mean   :224.5   Mean   :225  
#>  3rd Qu.:225.0   3rd Qu.:225  
#>  Max.   :225.0   Max.   :225  
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")

## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
#>       Asym       xmid      slope 
#> 258.440682   2.442061   1.858694 
fitted(mod1)
#>  [1]  94.34749 121.23271 137.45031 148.83053 157.45735 164.31866 169.95946
#>  [8] 174.71115 178.78954 182.34254 185.47566 188.26658 190.77402 193.04337
#> [15] 195.11033 197.00350 198.74606 200.35705 201.85227 203.24499 204.54643
#> [22] 205.76612 206.91229 207.99203 209.01150 209.97609 210.89054 211.75903
#> [29] 212.58527 213.37256 214.12386 214.84180 215.52877 216.18692 216.81820
#> [36] 217.42437 218.00703 218.56767 219.10762 219.62811 220.13027 220.61514
#> [43] 221.08369 221.53679 221.97528 222.39991 222.81138 223.21037 223.59747
#> [50] 223.97327
plot(sp1)
## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)

## Fit Arrhenius models to all random accumulations
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)

## Use nls() methods to the list of models
sapply(mods$models, AIC)
#>   [1] 339.6375 302.0388 311.1358 330.6704 326.7635 331.8803 310.5078 331.9992
#>   [9] 291.9614 316.1806 321.1436 330.3612 351.9384 366.4180 375.1539 314.8198
#>  [17] 365.6109 341.7665 336.8766 291.5889 313.0829 352.8956 365.4867 327.8955
#>  [25] 316.7024 334.7748 348.5523 309.3795 347.6068 361.8948 337.9602 346.4483
#>  [33] 345.8167 303.9609 327.4558 351.9462 306.6866 318.9341 336.6902 282.1333
#>  [41] 338.7179 315.2114 317.7183 345.3468 336.7774 313.0590 339.9524 344.0683
#>  [49] 341.7513 325.9581 330.4815 348.4942 308.9166 288.9574 336.6346 309.5134
#>  [57] 298.2728 343.7922 368.3111 322.5828 308.6839 315.4470 330.1068 346.9806
#>  [65] 337.2943 344.2560 364.8657 337.4751 330.5775 327.8712 340.3932 350.8871
#>  [73] 350.6292 339.5598 351.4364 296.5324 337.6837 334.7863 346.0071 374.7693
#>  [81] 373.7538 319.2654 343.6076 337.4138 330.2078 356.3056 301.1031 341.3953
#>  [89] 345.4966 342.5015 307.5547 334.1659 334.6058 308.5726 314.4447 370.6327
#>  [97] 311.9563 330.9361 312.4905 357.1784
```
