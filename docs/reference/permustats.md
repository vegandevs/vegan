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
#> Model    2.9966 5.4746 1.0147       0.9504 1.6593    0.001 ***
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
#> Overall (F)    1.9506  0.6068  1.2355          0.8981  2.6141    0.183  
#> BF-HF (t)     -0.5634 -0.4360 -0.0011 -2.0295 -0.0547  2.1039    0.607  
#> BF-NM (t)     -2.2387 -1.7794  0.0109 -2.1236  0.0337  2.0142    0.081 .
#> BF-SF (t)     -1.1675 -0.9180 -0.0303 -1.9206 -0.0247  1.7872    0.323  
#> HF-NM (t)     -2.1017 -1.8570  0.0562 -1.8001  0.0389  1.9507    0.059 .
#> HF-SF (t)     -0.8789 -0.7560  0.0041 -1.9673 -0.0375  1.8828    0.393  
#> NM-SF (t)      0.9485  0.8133 -0.0297 -1.9023 -0.0136  1.8619    0.383  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Interval (Upper - Lower) = 0.9)
```
