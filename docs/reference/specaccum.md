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
#> Sites     1.00000   2.00000   3.00000   4.0000   5.00000   6.00000   7.00000
#> Richness 91.20000 121.31000 137.79000 150.0700 158.94000 165.91000 171.32000
#> sd        7.56387   7.24701   7.18766   6.7678   6.44138   5.78765   5.65485
#>                                                                            
#> Sites      8.00000   9.00000  10.000  11.0000  12.00000  13.00000  14.00000
#> Richness 175.72000 179.81000 183.090 185.5200 188.07000 190.65000 192.71000
#> sd         5.35239   4.90453   4.816   4.7022   4.86018   4.78925   4.53114
#>                                                                              
#> Sites     15.00000  16.0000  17.00000  18.00000  19.00000  20.00000  21.00000
#> Richness 194.63000 196.3100 198.05000 199.63000 201.15000 202.49000 203.79000
#> sd         4.47362   4.2822   4.06109   3.93804   3.82806   3.55759   3.55986
#>                                                                              
#> Sites     22.0000  23.00000  24.00000  25.00000  26.00000  27.00000  28.00000
#> Richness 204.9400 206.23000 207.29000 208.42000 209.50000 210.58000 211.54000
#> sd         3.6703   3.60094   3.59094   3.60185   3.42746   3.43829   3.43899
#>                                                                               
#> Sites     29.00000  30.00000  31.00000  32.00000  33.00000  34.00000  35.00000
#> Richness 212.46000 213.25000 214.11000 214.87000 215.75000 216.42000 217.04000
#> sd         3.25179   3.06619   3.01141   2.97686   2.99284   2.79314   2.73001
#>                                                                               
#> Sites     36.00000  37.00000  38.00000  39.00000  40.00000  41.00000  42.00000
#> Richness 217.68000 218.23000 218.92000 219.56000 220.25000 220.87000 221.29000
#> sd         2.60101   2.52204   2.39393   2.28884   2.15732   2.15395   2.00653
#>                                                                               
#> Sites     43.00000  44.00000  45.00000  46.00000  47.00000  48.00000  49.00000
#> Richness 221.74000 222.28000 222.77000 223.12000 223.58000 224.17000 224.62000
#> sd         1.79573   1.66412   1.55606   1.35795   1.28062   0.89955   0.63214
#>             
#> Sites     50
#> Richness 225
#> sd         0
summary(sp2)
#>  1 sites          2 sites         3 sites         4 sites        
#>  Min.   : 80.00   Min.   :103.0   Min.   :122.0   Min.   :135.0  
#>  1st Qu.: 85.00   1st Qu.:116.0   1st Qu.:133.0   1st Qu.:145.0  
#>  Median : 90.00   Median :121.0   Median :138.0   Median :150.0  
#>  Mean   : 91.20   Mean   :121.3   Mean   :137.8   Mean   :150.1  
#>  3rd Qu.: 95.75   3rd Qu.:127.0   3rd Qu.:143.0   3rd Qu.:155.2  
#>  Max.   :109.00   Max.   :140.0   Max.   :153.0   Max.   :166.0  
#>  5 sites         6 sites         7 sites         8 sites        
#>  Min.   :142.0   Min.   :151.0   Min.   :156.0   Min.   :164.0  
#>  1st Qu.:154.0   1st Qu.:162.0   1st Qu.:168.0   1st Qu.:173.0  
#>  Median :158.5   Median :165.5   Median :171.5   Median :175.0  
#>  Mean   :158.9   Mean   :165.9   Mean   :171.3   Mean   :175.7  
#>  3rd Qu.:164.0   3rd Qu.:170.0   3rd Qu.:175.2   3rd Qu.:180.0  
#>  Max.   :173.0   Max.   :179.0   Max.   :185.0   Max.   :188.0  
#>  9 sites         10 sites        11 sites        12 sites       
#>  Min.   :168.0   Min.   :172.0   Min.   :174.0   Min.   :174.0  
#>  1st Qu.:176.0   1st Qu.:180.0   1st Qu.:182.0   1st Qu.:186.0  
#>  Median :180.0   Median :183.0   Median :185.0   Median :188.0  
#>  Mean   :179.8   Mean   :183.1   Mean   :185.5   Mean   :188.1  
#>  3rd Qu.:183.0   3rd Qu.:186.0   3rd Qu.:189.0   3rd Qu.:192.0  
#>  Max.   :194.0   Max.   :194.0   Max.   :195.0   Max.   :198.0  
#>  13 sites        14 sites        15 sites        16 sites       
#>  Min.   :177.0   Min.   :179.0   Min.   :182.0   Min.   :184.0  
#>  1st Qu.:188.0   1st Qu.:190.0   1st Qu.:191.0   1st Qu.:194.0  
#>  Median :191.0   Median :194.0   Median :195.0   Median :197.0  
#>  Mean   :190.7   Mean   :192.7   Mean   :194.6   Mean   :196.3  
#>  3rd Qu.:194.0   3rd Qu.:195.2   3rd Qu.:198.0   3rd Qu.:199.0  
#>  Max.   :201.0   Max.   :202.0   Max.   :204.0   Max.   :204.0  
#>  17 sites        18 sites        19 sites        20 sites       
#>  Min.   :185.0   Min.   :186.0   Min.   :187.0   Min.   :190.0  
#>  1st Qu.:195.8   1st Qu.:198.0   1st Qu.:199.0   1st Qu.:200.0  
#>  Median :198.0   Median :200.0   Median :201.0   Median :202.5  
#>  Mean   :198.1   Mean   :199.6   Mean   :201.2   Mean   :202.5  
#>  3rd Qu.:201.0   3rd Qu.:202.0   3rd Qu.:204.0   3rd Qu.:205.0  
#>  Max.   :207.0   Max.   :208.0   Max.   :210.0   Max.   :210.0  
#>  21 sites        22 sites        23 sites        24 sites       
#>  Min.   :192.0   Min.   :194.0   Min.   :195.0   Min.   :196.0  
#>  1st Qu.:202.0   1st Qu.:202.0   1st Qu.:204.0   1st Qu.:205.0  
#>  Median :204.0   Median :205.0   Median :207.0   Median :208.0  
#>  Mean   :203.8   Mean   :204.9   Mean   :206.2   Mean   :207.3  
#>  3rd Qu.:206.2   3rd Qu.:207.0   3rd Qu.:209.0   3rd Qu.:210.0  
#>  Max.   :211.0   Max.   :213.0   Max.   :213.0   Max.   :214.0  
#>  25 sites        26 sites        27 sites        28 sites       
#>  Min.   :198.0   Min.   :199.0   Min.   :200.0   Min.   :200.0  
#>  1st Qu.:206.0   1st Qu.:208.0   1st Qu.:209.0   1st Qu.:209.8  
#>  Median :209.0   Median :210.0   Median :210.5   Median :212.0  
#>  Mean   :208.4   Mean   :209.5   Mean   :210.6   Mean   :211.5  
#>  3rd Qu.:211.0   3rd Qu.:212.0   3rd Qu.:213.0   3rd Qu.:214.0  
#>  Max.   :217.0   Max.   :217.0   Max.   :218.0   Max.   :219.0  
#>  29 sites        30 sites        31 sites        32 sites       
#>  Min.   :201.0   Min.   :202.0   Min.   :203.0   Min.   :204.0  
#>  1st Qu.:211.0   1st Qu.:211.0   1st Qu.:213.0   1st Qu.:213.0  
#>  Median :213.0   Median :213.0   Median :214.0   Median :215.0  
#>  Mean   :212.5   Mean   :213.2   Mean   :214.1   Mean   :214.9  
#>  3rd Qu.:215.0   3rd Qu.:215.0   3rd Qu.:216.0   3rd Qu.:217.0  
#>  Max.   :219.0   Max.   :219.0   Max.   :220.0   Max.   :221.0  
#>  33 sites        34 sites        35 sites      36 sites        37 sites       
#>  Min.   :204.0   Min.   :208.0   Min.   :210   Min.   :211.0   Min.   :212.0  
#>  1st Qu.:214.0   1st Qu.:215.0   1st Qu.:216   1st Qu.:216.0   1st Qu.:217.0  
#>  Median :216.0   Median :216.0   Median :217   Median :218.0   Median :218.0  
#>  Mean   :215.8   Mean   :216.4   Mean   :217   Mean   :217.7   Mean   :218.2  
#>  3rd Qu.:218.0   3rd Qu.:218.2   3rd Qu.:219   3rd Qu.:220.0   3rd Qu.:220.0  
#>  Max.   :222.0   Max.   :223.0   Max.   :223   Max.   :224.0   Max.   :224.0  
#>  38 sites        39 sites        40 sites        41 sites       
#>  Min.   :213.0   Min.   :213.0   Min.   :214.0   Min.   :215.0  
#>  1st Qu.:217.0   1st Qu.:218.0   1st Qu.:219.0   1st Qu.:219.0  
#>  Median :219.0   Median :220.0   Median :220.0   Median :221.0  
#>  Mean   :218.9   Mean   :219.6   Mean   :220.2   Mean   :220.9  
#>  3rd Qu.:220.0   3rd Qu.:221.0   3rd Qu.:222.0   3rd Qu.:222.0  
#>  Max.   :225.0   Max.   :225.0   Max.   :225.0   Max.   :225.0  
#>  42 sites        43 sites        44 sites        45 sites       
#>  Min.   :216.0   Min.   :216.0   Min.   :217.0   Min.   :217.0  
#>  1st Qu.:220.0   1st Qu.:221.0   1st Qu.:221.0   1st Qu.:222.0  
#>  Median :221.5   Median :222.0   Median :222.0   Median :223.0  
#>  Mean   :221.3   Mean   :221.7   Mean   :222.3   Mean   :222.8  
#>  3rd Qu.:223.0   3rd Qu.:223.0   3rd Qu.:223.2   3rd Qu.:224.0  
#>  Max.   :225.0   Max.   :225.0   Max.   :225.0   Max.   :225.0  
#>  46 sites        47 sites        48 sites        49 sites        50 sites     
#>  Min.   :218.0   Min.   :219.0   Min.   :222.0   Min.   :222.0   Min.   :225  
#>  1st Qu.:222.0   1st Qu.:223.0   1st Qu.:224.0   1st Qu.:224.0   1st Qu.:225  
#>  Median :223.0   Median :224.0   Median :224.0   Median :225.0   Median :225  
#>  Mean   :223.1   Mean   :223.6   Mean   :224.2   Mean   :224.6   Mean   :225  
#>  3rd Qu.:224.0   3rd Qu.:224.0   3rd Qu.:225.0   3rd Qu.:225.0   3rd Qu.:225  
#>  Max.   :225.0   Max.   :225.0   Max.   :225.0   Max.   :225.0   Max.   :225  
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
#>   [1] 347.2816 325.6343 315.9136 345.4793 297.2704 342.9443 343.5358 308.3635
#>   [9] 345.7665 332.2003 327.9394 336.3936 345.3200 347.5309 313.0432 336.3863
#>  [17] 316.4203 317.1504 327.3126 332.6625 342.4624 351.3331 342.5683 336.9944
#>  [25] 321.5591 360.5222 320.7123 310.1364 345.6802 286.4162 313.5608 336.6158
#>  [33] 322.3123 296.6553 348.3919 335.0432 301.8475 362.5773 361.6093 326.3865
#>  [41] 335.7143 343.5142 359.0284 303.7928 329.4691 334.8476 312.2365 303.6824
#>  [49] 327.0915 352.2847 353.3203 315.8118 362.1672 326.9056 349.2616 342.0743
#>  [57] 345.1636 341.7187 315.1236 354.5164 332.7206 336.3892 338.0730 332.4325
#>  [65] 346.6708 315.8174 352.8257 286.1445 346.6942 329.1196 341.1443 314.8350
#>  [73] 336.5798 357.8063 328.1803 326.7973 329.7093 329.7821 318.1719 323.7645
#>  [81] 288.7501 361.3951 324.2064 319.1677 346.1219 349.4702 361.0926 326.3232
#>  [89] 342.1199 347.5353 302.2238 339.8217 333.2684 341.4940 334.9774 349.7622
#>  [97] 348.2272 351.9914 329.0865 342.1391
```
