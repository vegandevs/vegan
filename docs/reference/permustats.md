# Extract, Analyse and Display Permutation Results

The `permustats` function extracts permutation results of vegan
functions. Its support functions can find quantiles and standardized
effect sizes, plot densities and Q-Q plots.

## Usage

``` r
permustats(x, ...)
# S3 method for class 'permustats'
summary(object, interval = 0.95, alternative, ...)
# S3 method for class 'permustats'
density(x, observed = TRUE, ...)
# S3 method for class 'permustats'
qqnorm(y, observed = TRUE, ...)
# S3 method for class 'permustats'
boxplot(x, scale = FALSE, names, ...)
# S3 method for class 'permustats'
pairs(x, ...)
```

## Arguments

- object, x, y:

  The object to be handled.

- interval:

  numeric; the coverage interval reported.

- alternative:

  A character string specifying the limits used for the `interval` and
  the direction of the test when evaluating the \\p\\-values. Must be
  one of `"two.sided"` (both upper and lower limit), `"greater"` (upper
  limit), `"less"` (lower limit). Usually `alternative` is given in the
  result object, but it can be specified with this argument.

- observed:

  Add observed statistic among permutations.

- scale:

  Use standardized effect size (SES).

- names:

  Names of boxes (default: names of statistics).

- ...:

  Other arguments passed to the function. In `density` these are passed
  to [`density.default`](https://rdrr.io/r/stats/density.html), and in
  `boxplot` to
  [`boxplot.default`](https://rdrr.io/r/graphics/boxplot.html).

## Details

The `permustats` function extracts permutation results and observed
statistics from several vegan functions that perform permutations or
simulations.

The `summary` method of `permustats` estimates the standardized effect
sizes (SES) as the difference of observed statistic and mean of
permutations divided by the standard deviation of permutations (also
known as \\z\\-values). It also prints the the mean, median, and limits
which contain `interval` percent of permuted values. With the default
(`interval = 0.95`), for two-sided test these are (2.5%, 97.5%) and for
one-sided tests either 5% or 95% quantile and the \\p\\-value depending
on the test direction. The mean, quantiles and \\z\\ values are
evaluated from permuted values without observed statistic, but the
\\p\\-value is evaluated with the observed statistic. The intervals and
the \\p\\-value are evaluated with the same test direction as in the
original test, but this can be changed with argument `alternative`.
Several `permustats` objects can be combined with `c` function. The `c`
function checks that statistics are equal, but performs no other sanity
tests.

The results can be displayed with conventional graphics or as ggplot2
graphics using `autoplot` function in
[ggvegan](https://CRAN.R-project.org/package=ggvegan) package.

The `density` and `densityplot` methods display the kernel density
estimates of permuted values. When observed value of the statistic is
included in the permuted values, the `densityplot` method marks the
observed statistic as a vertical line. However the `density` method uses
its standard `plot` method and cannot mark the observed value. Only one
statistic can be displayed with `density` and for several statistics
`permulattice` or `densityplot` must be used.

The `qqnorm` method display Q-Q plots of permutations, optionally
together with the observed value (default) which is shown as horizontal
line in plots. `qqnorm` plots permutation values against standard Normal
variate. The permutations are standardized without the observed
statistic, similarly as in `summary`. Only one statistic can be shown
with `qqnorm` and for several statistics ggvegan package must be used.

Function `boxplot` draws the box-and-whiskers plots of effect size, or
the difference of permutations and observed statistic. If
`scale = TRUE`, permutations are standardized to unit standard
deviation, and the plot will show the standardized effect sizes.

Function `pairs` plots permutation values of statistics against each
other. The function passes extra arguments to
[`pairs`](https://rdrr.io/r/graphics/pairs.html).

The `permustats` can extract permutation statistics from the results of
[`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.md),
[`anosim`](https://vegandevs.github.io/vegan/reference/anosim.md),
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`mantel`](https://vegandevs.github.io/vegan/reference/mantel.md),
[`mantel.partial`](https://vegandevs.github.io/vegan/reference/mantel.md),
[`mrpp`](https://vegandevs.github.io/vegan/reference/mrpp.md),
[`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md),
[`ordiareatest`](https://vegandevs.github.io/vegan/reference/ordihull.md),
[`permutest.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`protest`](https://vegandevs.github.io/vegan/reference/procrustes.md),
and
[`permutest.betadisper`](https://vegandevs.github.io/vegan/reference/permutest.betadisper.md).

## Value

The `permustats` function returns an object of class `"permustats"`.
This is a list of items `"statistic"` for observed statistics,
`permutations` which contains permuted values, and `alternative` which
contains text defining the character of the test (`"two.sided"`,
`"less"` or `"greater"`). The
[`qqnorm`](https://rdrr.io/r/stats/qqnorm.html) and
[`density`](https://rdrr.io/r/stats/density.html) methods return their
standard result objects.

## Author

Jari Oksanen with contributions from Gavin L. Simpson
(`permustats.permutest.betadisper` method and related modifications to
`summary.permustats` and the `print` method) and Eduard Szöcs
(`permustats.anova.cca).`

## See also

[`density`](https://rdrr.io/r/stats/density.html),
[`qqnorm`](https://rdrr.io/r/stats/qqnorm.html),
[`boxplot`](https://rdrr.io/r/graphics/boxplot.html).

## Examples

``` r
data(dune, dune.env)
mod <- adonis2(dune ~ Management + A1, data = dune.env)
## use permustats
perm <- permustats(mod)
summary(perm)
#> 
#>       statistic    SES   mean lower median  upper Pr(perm)    
#> Model    2.9966 5.4237 1.0137       0.9495 1.6706    0.001 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Interval (Upper - Lower) = 0.95)
boxplot(perm, scale=TRUE, lty=1, pch=16, cex=0.6, col="hotpink", ylab="SES")
abline(h=0, col="skyblue")

## example of multiple types of statistic
mod <- with(dune.env, betadisper(vegdist(dune), Management))
pmod <- permutest(mod, nperm = 99, pairwise = TRUE)
perm <- permustats(pmod)
summary(perm, interval = 0.90)
#> 
#>             statistic     SES    mean   lower  median   upper Pr(perm)  
#> Overall (F)    1.9506  0.6077  1.2335          0.8982  2.5872    0.183  
#> BF-HF (t)     -0.5634 -0.4466  0.0187 -2.0376 -0.0157  2.1457    0.605  
#> BF-NM (t)     -2.2387 -1.7986  0.0197 -2.0642  0.0297  2.0142    0.073 .
#> BF-SF (t)     -1.1675 -0.9350 -0.0249 -1.9101 -0.0277  1.7764    0.315  
#> HF-NM (t)     -2.1017 -1.8283  0.0365 -1.8312  0.0279  1.9724    0.063 .
#> HF-SF (t)     -0.8789 -0.7373 -0.0144 -1.9970 -0.0642  1.8540    0.405  
#> NM-SF (t)      0.9485  0.8218 -0.0305 -1.9023 -0.0214  1.8619    0.371  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Interval (Upper - Lower) = 0.9)
```
