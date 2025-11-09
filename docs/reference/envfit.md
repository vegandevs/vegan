# Fits an Environmental Vector or Factor onto an Ordination

The function fits environmental vectors or factors onto an ordination.
The projections of points onto vectors have maximum correlation with
corresponding environmental variables, and the factors show the averages
of factor levels. For continuous varaibles this is equal to fitting a
linear trend surface (plane in 2D) for a variable (see
[`ordisurf`](https://vegandevs.github.io/vegan/reference/ordisurf.md));
this trend surface can be presented by showing its gradient (direction
of steepest increase) using an arrow. The environmental variables are
the dependent variables that are explained by the ordination scores, and
each dependent variable is analysed separately.

## Usage

``` r
# Default S3 method
envfit(ord, env, permutations = 999, strata = NULL, 
   choices=c(1,2),  display = "sites", w, na.rm = FALSE, ...)
# S3 method for class 'formula'
envfit(formula, data, ...)
# S3 method for class 'envfit'
plot(x, choices = c(1,2), labels, arrow.mul, at = c(0,0), 
   axis = FALSE, p.max = NULL, r2.min = NULL, col = "blue", bg, add = TRUE, ...)
# S3 method for class 'envfit'
scores(x, display, choices, arrow.mul=1, tidy = FALSE, ...)
vectorfit(X, P, permutations = 0, strata = NULL, w, ...)
factorfit(X, P, permutations = 0, strata = NULL, w, ...)
```

## Arguments

- ord:

  An ordination object or other structure from which the ordination
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md) can
  be extracted (including a data frame or matrix of scores).

- env:

  Data frame, matrix or vector of environmental variables. The variables
  can be of mixed type (factors, continuous variables) in data frames.

- X:

  Matrix or data frame of ordination scores.

- P:

  Data frame, matrix or vector of environmental variable(s). These must
  be continuous for `vectorfit` and factors or characters for
  `factorfit`.

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required, or a permutation matrix where each
  row gives the permuted indices. Set `permutations = 0` to skip
  permutations.

- formula, data:

  Model [`formula`](https://rdrr.io/r/stats/formula.html) and data.

- na.rm:

  Remove points with missing values in ordination scores or
  environmental variables. The operation is casewise: the whole row of
  data is removed if there is a missing value and `na.rm = TRUE`.

- x:

  A result object from `envfit`. For `ordiArrowMul` and
  `ordiArrowTextXY` this must be a two-column matrix (or matrix-like
  object) containing the coordinates of arrow heads on the two plot
  axes, and other methods extract such a structure from the `envfit`
  results.

- choices:

  Axes to plotted.

- tidy:

  Return scores that are compatible with
  [ggplot2](https://CRAN.R-project.org/package=ggplot2): all scores are
  in a single `data.frame`, score type is identified by factor variable
  `scores` (`"vectors"` or `"factors"`), the names by variable `label`.
  These scores are incompatible with conventional `plot` functions, but
  they can be used in ggplot2.

- labels:

  Change plotting labels. The argument should be a list with elements
  `vectors` and `factors` which give the new plotting labels. If either
  of these elements is omitted, the default labels will be used. If
  there is only one type of elements (only `vectors` or only `factors`),
  the labels can be given as vector. The default labels can be displayed
  with `labels` command.

- arrow.mul:

  Multiplier for vector lengths. The arrows are automatically scaled
  similarly as in
  [`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
  if this is not given in `plot` and `add = TRUE`. However, in `scores`
  it can be used to adjust arrow lengths when the `plot` function is not
  used.

- at:

  The origin of fitted arrows in the plot. If you plot arrows in other
  places then origin, you probably have to specify `arrrow.mul`.

- axis:

  Plot axis showing the scaling of fitted arrows.

- p.max, r2.min:

  Maximum estimated \\P\\ value and minimum \\r^2\\ for displayed
  variables. You must calculate \\P\\ values with setting `permutations`
  to use `p.max`.

- col:

  Colour in plotting.

- bg:

  Background colour for labels. If `bg` is set, the labels are displayed
  with
  [`ordilabel`](https://vegandevs.github.io/vegan/reference/ordilabel.md)
  instead of `text`. See Examples for using semitransparent background.

- add:

  Results added to an existing ordination plot.

- strata:

  An integer vector or factor specifying the strata for permutation. If
  supplied, observations are permuted only within the specified strata.

- display:

  In fitting functions these are ordinary site scores or linear
  combination scores (`"lc"`) in constrained ordination
  ([`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md)). In
  `scores` function they are either `"vectors"` or `"factors"` (with
  synonyms `"bp"` or `"cn"`, resp.).

- w:

  Weights used in fitting (concerns mainly
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md) and
  [`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)
  results which have nonconstant weights).

- ...:

  Parameters passed to
  [`scores`](https://vegandevs.github.io/vegan/reference/scores.md).

## Details

Function `envfit` finds vectors or factor averages of environmental
variables. Function `plot.envfit` adds these in an ordination diagram.
If `X` is a [`data.frame`](https://rdrr.io/r/base/data.frame.html),
`envfit` uses `factorfit` for
[`factor`](https://rdrr.io/r/base/factor.html) variables and `vectorfit`
for other variables. If `X` is a matrix or a vector, `envfit` uses only
`vectorfit`. Alternatively, the model can be defined a simplified model
[`formula`](https://rdrr.io/r/stats/formula.html), where the left hand
side must be an ordination result object or a matrix of ordination
scores, and right hand side lists the environmental variables. The
formula interface can be used for easier selection and/or transformation
of environmental variables. Only the main effects will be analysed even
if interaction terms were defined in the formula.

The ordination results are extracted with
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md) and
all extra arguments are passed to the `scores`. The fitted models only
apply to the results defined when extracting the scores when using
`envfit`. For instance, `scaling` in constrained ordination (see
[`scores.rda`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md))
must be set in the same way in `envfit` and in the `plot` or the
ordination results (see Examples).

The printed output of continuous variables (vectors) gives the direction
cosines which are the coordinates of the heads of unit length vectors.
In `plot` these are scaled by their correlation (square root of the
column `r2`) so that “weak” predictors have shorter arrows than “strong”
predictors. You can see the scaled relative lengths using command
`scores`. The `plot`ted (and scaled) arrows are further adjusted to the
current graph using a constant multiplier: this will keep the relative
`r2`-scaled lengths of the arrows but tries to fill the current plot.
You can see the multiplier using `ordiArrowMul(result_of_envfit)`, and
set it with the argument `arrow.mul`.

Functions `vectorfit` and `factorfit` can be called directly. Function
`vectorfit` finds directions in the ordination space towards which the
environmental vectors change most rapidly and to which they have maximal
correlations with the ordination configuration. Function `factorfit`
finds averages of ordination scores for factor levels. Function
`factorfit` treats ordered and unordered factors similarly.

If `permutations` \\\> 0\\, the significance of fitted vectors or
factors is assessed using permutation of environmental variables. The
goodness of fit statistic is squared correlation coefficient (\\r^2\\).
For factors this is defined as \\r^2 = 1 - ss_w/ss_t\\, where \\ss_w\\
and \\ss_t\\ are within-group and total sums of squares. See
[`permutations`](https://vegandevs.github.io/vegan/reference/permutations.md)
for additional details on permutation tests in Vegan.

User can supply a vector of prior weights `w`. If the ordination object
has weights, these will be used. In practise this means that the row
totals are used as weights with
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) or
[`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)
results. If you do not like this, but want to give equal weights to all
sites, you should set `w = NULL`. The fitted vectors are similar to
biplot arrows in constrained ordination only when fitted to LC scores
(`display = "lc"`) and you set `scaling = "species"` (see
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)).
The weighted fitting gives similar results to biplot arrows and class
centroids in
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md).

The lengths of arrows for fitted vectors are automatically adjusted for
the physical size of the plot, and the arrow lengths cannot be compared
across plots. For similar scaling of arrows, you must explicitly set the
`arrow.mul` argument in the `plot` command; see
[`ordiArrowMul`](https://vegandevs.github.io/vegan/reference/ordiArrowTextXY.md)
and
[`ordiArrowTextXY`](https://vegandevs.github.io/vegan/reference/ordiArrowTextXY.md).

The results can be accessed with `scores.envfit` function which returns
either the fitted vectors scaled by correlation coefficient or the
centroids of the fitted environmental variables, or a named list of
both.

## Value

Functions `vectorfit` and `factorfit` return lists of classes
`vectorfit` and `factorfit` which have a `print` method. The result
object have the following items:

- arrows:

  Arrow endpoints from `vectorfit`. The arrows are scaled to unit
  length.

- centroids:

  Class centroids from `factorfit`.

- r:

  Goodness of fit statistic: Squared correlation coefficient

- permutations:

  Number of permutations.

- control:

  A list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html).

- pvals:

  Empirical P-values for each variable.

Function `envfit` returns a list of class `envfit` with results of
`vectorfit` and `envfit` as items.

Function `plot.envfit` scales the vectors by correlation.

## Author

Jari Oksanen. The permutation test derives from the code suggested by
Michael Scroggie.

## Note

Fitted vectors have become the method of choice in displaying
environmental variables in ordination. Indeed, they are the optimal way
of presenting environmental variables in Constrained Correspondence
Analysis [`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
since there they are the linear constraints. In unconstrained ordination
the relation between external variables and ordination configuration may
be less linear, and therefore other methods than arrows may be more
useful. The simplest is to adjust the plotting symbol sizes (`cex`,
[`symbols`](https://rdrr.io/r/graphics/symbols.html)) by environmental
variables. Fancier methods involve smoothing and regression methods that
abound in R, and
[`ordisurf`](https://vegandevs.github.io/vegan/reference/ordisurf.md)
provides a wrapper for some.

## See also

A better alternative to vectors may be
[`ordisurf`](https://vegandevs.github.io/vegan/reference/ordisurf.md).

## Examples

``` r
data(varespec, varechem)
library(MASS)
ord <- metaMDS(varespec)
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.1843196 
#> Run 1 stress 0.1955836 
#> Run 2 stress 0.18458 
#> ... Procrustes: rmse 0.0493631  max resid 0.1575579 
#> Run 3 stress 0.1976156 
#> Run 4 stress 0.2094754 
#> Run 5 stress 0.1974418 
#> Run 6 stress 0.2185648 
#> Run 7 stress 0.2092456 
#> Run 8 stress 0.18584 
#> Run 9 stress 0.2394817 
#> Run 10 stress 0.2087932 
#> Run 11 stress 0.1843196 
#> ... New best solution
#> ... Procrustes: rmse 2.517178e-05  max resid 0.0001011706 
#> ... Similar to previous best
#> Run 12 stress 0.2240541 
#> Run 13 stress 0.1974419 
#> Run 14 stress 0.195049 
#> Run 15 stress 0.1962451 
#> Run 16 stress 0.1869637 
#> Run 17 stress 0.1843196 
#> ... Procrustes: rmse 2.289873e-05  max resid 8.006074e-05 
#> ... Similar to previous best
#> Run 18 stress 0.1845801 
#> ... Procrustes: rmse 0.04937852  max resid 0.15765 
#> Run 19 stress 0.1967393 
#> Run 20 stress 0.2028828 
#> *** Best solution repeated 2 times
(fit <- envfit(ord, varechem, perm = 999))
#> 
#> ***VECTORS
#> 
#>             NMDS1    NMDS2     r2 Pr(>r)    
#> N        -0.05039 -0.99873 0.2081  0.080 .  
#> P         0.68701  0.72665 0.1755  0.147    
#> K         0.82728  0.56179 0.1657  0.171    
#> Ca        0.75014  0.66128 0.2811  0.047 *  
#> Mg        0.69675  0.71732 0.3494  0.021 *  
#> S         0.27625  0.96109 0.1774  0.133    
#> Al       -0.83778  0.54601 0.5155  0.001 ***
#> Fe       -0.86199  0.50692 0.4001  0.003 ** 
#> Mn        0.80234 -0.59686 0.5322  0.002 ** 
#> Zn        0.66517  0.74669 0.1779  0.142    
#> Mo       -0.84876  0.52877 0.0517  0.572    
#> Baresoil  0.87211 -0.48931 0.2494  0.052 .  
#> Humdepth  0.92637 -0.37662 0.5590  0.001 ***
#> pH       -0.79908  0.60123 0.2624  0.052 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Permutation: free
#> Number of permutations: 999
#> 
#> 
scores(fit, "vectors")
#>                NMDS1      NMDS2
#> N        -0.02298558 -0.4555563
#> P         0.28784302  0.3044503
#> K         0.33676522  0.2286909
#> Ca        0.39771862  0.3506023
#> Mg        0.41183998  0.4239993
#> S         0.11635223  0.4047925
#> Al       -0.60152079  0.3920315
#> Fe       -0.54522417  0.3206335
#> Mn        0.58534888 -0.4354406
#> Zn        0.28057323  0.3149614
#> Mo       -0.19295095  0.1202064
#> Baresoil  0.43555136 -0.2443704
#> Humdepth  0.69263050 -0.2815923
#> pH       -0.40930953  0.3079660
plot(ord)
plot(fit)
plot(fit, p.max = 0.05, col = "red")

## Adding fitted arrows to CCA. We use "lc" scores, and hope
## that arrows are scaled similarly in cca and envfit plots
ord <- cca(varespec ~ Al + P + K, varechem)
plot(ord, type="p")
fit <- envfit(ord, varechem, perm = 999, display = "lc")
plot(fit, p.max = 0.05, col = "red")

## 'scaling' must be set similarly in envfit and in ordination plot
plot(ord, type = "p", scaling = "sites")
fit <- envfit(ord, varechem, perm = 0, display = "lc", scaling = "sites")
plot(fit, col = "red")


## Class variables, formula interface, and displaying the
## inter-class variability with ordispider, and semitransparent
## white background for labels (semitransparent colours are not
## supported by all graphics devices)
data(dune)
data(dune.env)
ord <- cca(dune)
fit <- envfit(ord ~ Moisture + A1, dune.env, perm = 0)
plot(ord, type = "n")
with(dune.env, ordispider(ord, Moisture, col="skyblue"))
with(dune.env, points(ord, display = "sites", col = as.numeric(Moisture),
                      pch=16))
plot(fit, cex=1.2, axis=TRUE, bg = rgb(1, 1, 1, 0.5))

## Use shorter labels for factor centroids
labels(fit)
#> $vectors
#> [1] "A1"
#> 
#> $factors
#> [1] "Moisture1" "Moisture2" "Moisture4" "Moisture5"
#> 
plot(ord)
plot(fit, labels=list(factors = paste("M", c(1,2,4,5), sep = "")),
   bg = rgb(1,1,0,0.5))
```
