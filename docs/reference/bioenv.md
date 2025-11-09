# Best Subset of Environmental Variables with Maximum (Rank) Correlation with Community Dissimilarities

Function finds the best subset of environmental variables, so that the
Euclidean distances of scaled environmental variables have the maximum
(rank) correlation with community dissimilarities.

## Usage

``` r
# Default S3 method
bioenv(comm, env, method = "spearman", index = "bray",
       upto = ncol(env), trace = FALSE, partial = NULL, 
       metric = c("euclidean", "mahalanobis", "manhattan", "gower"),
       parallel = getOption("mc.cores"), ...)
# S3 method for class 'formula'
bioenv(formula, data, ...)
bioenvdist(x, which = "best")
```

## Arguments

- comm:

  Community data frame or a dissimilarity object or a square matrix that
  can be interpreted as dissimilarities.

- env:

  Data frame of continuous environmental variables.

- method:

  The correlation method used in
  [`cor`](https://rdrr.io/r/stats/cor.html).

- index:

  The dissimilarity index used for community data (`comm`) in
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md).
  This is ignored if `comm` are dissimilarities.

- upto:

  Maximum number of parameters in studied subsets.

- formula, data:

  Model [`formula`](https://rdrr.io/r/stats/formula.html) and data.

- trace:

  Trace the calculations

- partial:

  Dissimilarities partialled out when inspecting variables in `env`.

- metric:

  Metric used for distances of environmental distances. See Details.

- parallel:

  Number of parallel processes or a predefined socket cluster. With
  `parallel = 1` uses ordinary, non-parallel processing. The parallel
  processing is done with parallel package.

- x:

  `bioenv` result object.

- which:

  The number of the model for which the environmental distances are
  evaluated, or the `"best"` model.

- ...:

  Other arguments passed to
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md).

## Details

The function calculates a community dissimilarity matrix using
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md).
Then it selects all possible subsets of environmental variables,
[`scale`](https://rdrr.io/r/base/scale.html)s the variables, and
calculates Euclidean distances for this subset using
[`dist`](https://rdrr.io/r/stats/dist.html). The function finds the
correlation between community dissimilarities and environmental
distances, and for each size of subsets, saves the best result. There
are \\2^p-1\\ subsets of \\p\\ variables, and an exhaustive search may
take a very, very, very long time (parameter `upto` offers a partial
relief).

The argument `metric` defines distances in the given set of
environmental variables. With `metric = "euclidean"`, the variables are
scaled to unit variance and Euclidean distances are calculated. With
`metric = "mahalanobis"`, the Mahalanobis distances are calculated: in
addition to scaling to unit variance, the matrix of the current set of
environmental variables is also made orthogonal (uncorrelated). With
`metric = "manhattan"`, the variables are scaled to unit range and
Manhattan distances are calculated, so that the distances are sums of
differences of environmental variables. With `metric = "gower"`, the
Gower distances are calculated using function
[`daisy`](https://rdrr.io/pkg/cluster/man/daisy.html). This allows also
using factor variables, but with continuous variables the results are
equal to `metric = "manhattan"`.

The function can be called with a model
[`formula`](https://rdrr.io/r/stats/formula.html) where the LHS is the
data matrix and RHS lists the environmental variables. The formula
interface is practical in selecting or transforming environmental
variables.

With argument `partial` you can perform “partial” analysis. The
partializing item must be a dissimilarity object of class
[`dist`](https://rdrr.io/r/stats/dist.html). The `partial` item can be
used with any correlation `method`, but it is strictly correct only for
Pearson.

Function `bioenvdist` recalculates the environmental distances used
within the function. The default is to calculate distances for the best
model, but the number of any model can be given.

Clarke & Ainsworth (1993) suggested this method to be used for selecting
the best subset of environmental variables in interpreting results of
nonmetric multidimensional scaling (NMDS). They recommended a parallel
display of NMDS of community dissimilarities and NMDS of Euclidean
distances from the best subset of scaled environmental variables. They
warned against the use of Procrustes analysis, but to me this looks like
a good way of comparing these two ordinations.

Clarke & Ainsworth wrote a computer program BIO-ENV giving the name to
the current function. Presumably BIO-ENV was later incorporated in
Clarke's PRIMER software (available for Windows). In addition, Clarke &
Ainsworth suggested a novel method of rank correlation which is not
available in the current function.

## Value

The function returns an object of class `bioenv` with a `summary`
method.

## References

Clarke, K. R & Ainsworth, M. 1993. A method of linking multivariate
community structure to environmental variables. *Marine Ecology Progress
Series*, 92, 205–219.

## Author

Jari Oksanen

## Note

If you want to study the ‘significance’ of `bioenv` results, you can use
function
[`mantel`](https://vegandevs.github.io/vegan/reference/mantel.md) or
[`mantel.partial`](https://vegandevs.github.io/vegan/reference/mantel.md)
which use the same definition of correlation. However, `bioenv`
standardizes environmental variables depending on the used metric, and
you must do the same in
[`mantel`](https://vegandevs.github.io/vegan/reference/mantel.md) for
comparable results (the standardized data are returned as item `x` in
the result object). It is safest to use `bioenvdist` to extract the
environmental distances that really were used within `bioenv`. NB.,
`bioenv` selects variables to maximize the Mantel correlation, and
significance tests based on *a priori* selection of variables are
biased.

## See also

[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md),
[`dist`](https://rdrr.io/r/stats/dist.html),
[`cor`](https://rdrr.io/r/stats/cor.html) for underlying routines,
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md) and
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md) for
ordination,
[`procrustes`](https://vegandevs.github.io/vegan/reference/procrustes.md)
for Procrustes analysis,
[`protest`](https://vegandevs.github.io/vegan/reference/procrustes.md)
for an alternative, and
[`rankindex`](https://vegandevs.github.io/vegan/reference/rankindex.md)
for studying alternatives to the default Bray-Curtis index.

## Examples

``` r
# The method is very slow for large number of possible subsets.
# Therefore only 6 variables in this example.
data(varespec)
data(varechem)
sol <- bioenv(wisconsin(varespec) ~ log(N) + P + K + Ca + pH + Al, varechem)
sol
#> 
#> Call:
#> bioenv(formula = wisconsin(varespec) ~ log(N) + P + K + Ca +      pH + Al, data = varechem) 
#> 
#> Subset of environmental variables with best correlation to community data.
#> 
#> Correlations:    spearman 
#> Dissimilarities: bray 
#> Metric:          euclidean 
#> 
#> Best model has 3 parameters (max. 6 allowed):
#> P Ca Al
#> with correlation  0.4004806 
#> 
## IGNORE_RDIFF_BEGIN
summary(sol)
#>                     size correlation
#> P                      1      0.2513
#> P Al                   2      0.4004
#> P Ca Al                3      0.4005
#> P Ca pH Al             4      0.3619
#> log(N) P Ca pH Al      5      0.3216
#> log(N) P K Ca pH Al    6      0.2822
## IGNORE_RDIFF_END
```
