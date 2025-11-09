# Principal Coordinates Analysis and \[Partial\] Distance-based Redundancy Analysis

Distance-based redundancy analysis (dbRDA) is an ordination method
similar to Redundancy Analysis
([`rda`](https://vegandevs.github.io/vegan/reference/cca.md)), but it
allows non-Euclidean dissimilarity indices, such as Manhattan or
Bray-Curtis distance. Despite this non-Euclidean feature, the analysis
is strictly linear and metric. If called with Euclidean distance, the
results are identical to
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md), but dbRDA
will be less efficient. Functions `dbrda` is constrained versions of
metric scaling, a.k.a. principal coordinates analysis, which are based
on the Euclidean distance but can be used, and are more useful, with
other dissimilarity measures. Function `capscale` is a simplified
version based on Euclidean approximation of dissimilarities. The
functions can also perform unconstrained principal coordinates analysis
(PCO), optionally using extended dissimilarities. `pco()` is a wrapper
to `dbrda()`, which performs PCO.

## Usage

``` r
dbrda(formula, data, distance = "euclidean", sqrt.dist = FALSE,
    add = FALSE, dfun = vegdist, metaMDSdist = FALSE,
    na.action = na.fail, subset = NULL, ...)
capscale(formula, data, distance = "euclidean", sqrt.dist = FALSE,
    comm = NULL, add = FALSE,  dfun = vegdist, metaMDSdist = FALSE,
    na.action = na.fail, subset = NULL, ...)
pco(X, ...)
```

## Arguments

- formula:

  Model formula. The function can be called only with the formula
  interface. Most usual features of
  [`formula`](https://rdrr.io/r/stats/formula.html) hold, especially as
  defined in [`cca`](https://vegandevs.github.io/vegan/reference/cca.md)
  and [`rda`](https://vegandevs.github.io/vegan/reference/cca.md). The
  LHS must be either a community data matrix or a dissimilarity matrix,
  e.g., from
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md) or
  [`dist`](https://rdrr.io/r/stats/dist.html). If the LHS is a data
  matrix, function
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md) or
  function given in `dfun` will be used to find the dissimilarities. The
  RHS defines the constraints. The constraints can be continuous
  variables or factors, they can be transformed within the formula, and
  they can have interactions as in a typical
  [`formula`](https://rdrr.io/r/stats/formula.html). The RHS can have a
  special term `Condition` that defines variables to be “partialled out”
  before constraints, just like in
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md) or
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md). This
  allows the use of partial dbRDA.

- X:

  Community data matrix.

- data:

  Data frame containing the variables on the right hand side of the
  model formula.

- distance:

  The name of the dissimilarity (or distance) index if the LHS of the
  `formula` is a data frame instead of dissimilarity matrix.

- sqrt.dist:

  Take square roots of dissimilarities. See section `Details` below.

- comm:

  Community data frame which will be used for finding species scores
  when the LHS of the `formula` was a dissimilarity matrix. This is not
  used if the LHS is a data frame. If this is not supplied, the “species
  scores” are unavailable when dissimilarities were supplied. N.B., this
  is only available in `capscale`: `dbrda` does not return species
  scores. Function
  [`sppscores`](https://vegandevs.github.io/vegan/reference/sppscores.md)
  can be used to add species scores if they are missing.

- add:

  Add a constant to the non-diagonal dissimilarities such that all
  eigenvalues are non-negative in the underlying Principal Co-ordinates
  Analysis (see
  [`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md)
  for details). `"lingoes"` (or `TRUE`) uses the recommended method of
  Legendre & Anderson (1999: “method 1”) and `"cailliez"` uses their
  “method 2”. The latter is the only one in
  [`cmdscale`](https://rdrr.io/r/stats/cmdscale.html).

- dfun:

  Distance or dissimilarity function used. Any function returning
  standard `"dist"` and taking the index name as the first argument can
  be used.

- metaMDSdist:

  Use
  [`metaMDSdist`](https://vegandevs.github.io/vegan/reference/metaMDS.md)
  similarly as in
  [`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md).
  This means automatic data transformation and using extended flexible
  shortest path dissimilarities (function
  [`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md))
  when there are many dissimilarities based on no shared species.

- na.action:

  Handling of missing values in constraints or conditions. The default
  ([`na.fail`](https://rdrr.io/r/stats/na.fail.html)) is to stop with
  missing values. Choices
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html) and
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) delete rows with
  missing values, but differ in representation of results. With
  `na.omit` only non-missing site scores are shown, but `na.exclude`
  gives `NA` for scores of missing observations. Unlike in
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md), no WA
  scores are available for missing constraints or conditions.

- subset:

  Subset of data rows. This can be a logical vector which is `TRUE` for
  kept observations, or a logical expression which can contain variables
  in the working environment, `data` or species names of the community
  data (if given in the formula or as `comm` argument).

- ...:

  Other parameters passed to underlying functions (e.g.,
  [`metaMDSdist`](https://vegandevs.github.io/vegan/reference/metaMDS.md)).
  For `pco()` argument are passed to `dbrda()`.

## Details

Functions `dbrda` and `capscale` provide two alternative implementations
of dbRDA. Function `dbrda` is based on McArdle & Anderson (2001) and
directly decomposes dissimilarities. With Euclidean distances results
are identical to
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md).
Non-Euclidean dissimilarities may give negative eigenvalues associated
with imaginary axes. Function `capscale` is based on Legendre & Anderson
(1999): the dissimilarity data are first ordinated using metric scaling,
and the ordination results are analysed as
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md). `capscale`
ignores the imaginary component and will not give negative eigenvalues
(but will report the magnitude on imaginary component).

If the user supplied a community data frame instead of dissimilarities,
the functions will find dissimilarities using
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md) or
distance function given in `dfun` with specified `distance`. The
functions will accept distance objects from
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md),
[`dist`](https://rdrr.io/r/stats/dist.html), or any other method
producing compatible objects. The constraining variables can be
continuous or factors or both, they can have interaction terms, or they
can be transformed in the call. Moreover, there can be a special term
`Condition` just like in
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) and
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) so that
“partial” analysis can be performed.

Function `dbrda` does not return species scores, and they can also be
missing in `capscale`, but they can be added after the analysis using
function
[`sppscores`](https://vegandevs.github.io/vegan/reference/sppscores.md).

Non-Euclidean dissimilarities can produce negative eigenvalues (Legendre
& Anderson 1999, McArdle & Anderson 2001). If there are negative
eigenvalues, the printed output of `capscale` will add a column with
sums of positive eigenvalues and an item of sum of negative eigenvalues,
and `dbrda` will add a column giving the number of real dimensions with
positive eigenvalues. If negative eigenvalues are disturbing, functions
let you distort the dissimilarities so that only non-negative
eigenvalues will be produced with argument `add = TRUE`. Alternatively,
with `sqrt.dist = TRUE`, square roots of dissimilarities can be used
which may help in avoiding negative eigenvalues (Legendre & Anderson
1999).

The functions can be also used to perform ordinary metric scaling a.k.a.
principal coordinates analysis by using a formula with only a constant
on the right hand side, or `comm ~ 1`. The new function `pco()`
implements principal coordinates analysis via `dbrda()` directly, using
this formula. With `metaMDSdist = TRUE`, the function can do automatic
data standardization and use extended dissimilarities using function
[`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md)
similarly as in non-metric multidimensional scaling with
[`metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.md).

## Value

The functions return an object of class `dbrda` or `capscale` which
inherit from
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md). See
[`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.md)
for description of the result object. Function `pco()` returns an object
of class `"vegan_pco"` (which inherits from class `"dbrda"`) to avoid
clashes with other packages.

## References

Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal
coordinates: a useful method of constrained ordination for ecology.
*Ecology* 84, 511–525.

Gower, J.C. (1985). Properties of Euclidean and non-Euclidean distance
matrices. *Linear Algebra and its Applications* 67, 81–97.

Legendre, P. & Anderson, M. J. (1999). Distance-based redundancy
analysis: testing multispecies responses in multifactorial ecological
experiments. *Ecological Monographs* 69, 1–24.

Legendre, P. & Legendre, L. (2012). *Numerical Ecology*. 3rd English
Edition. Elsevier.

McArdle, B.H. & Anderson, M.J. (2001). Fitting multivariate models to
community data: a comment on distance-based redundancy analysis.
*Ecology* 82, 290–297.

## Author

Jari Oksanen

## Note

Function `dbrda` implements real distance-based RDA and is preferred
over `capscale`. `capscale` was originally developed as a variant of
constrained analysis of proximities (Anderson & Willis 2003), but these
developments made it more similar to dbRDA. However, it discards the
imaginary dimensions with negative eigenvalues and ordination and
significance tests area only based on real dimensions and positive
eigenvalues. `capscale` may be removed from vegan in the future. It has
been in `vegan` since 2003 (CRAN release 1.6-0) while `dbrda` was first
released in 2016 (version 2.4-0), and removal of `capscale` may be
disruptive to historical examples and scripts, but in modern times
`dbrda` should be used.

The inertia is named after the dissimilarity index as defined in the
dissimilarity data, or as `unknown distance` if such information is
missing. If the largest original dissimilarity was larger than 4,
`capscale` handles input similarly as `rda` and bases its analysis on
variance instead of sum of squares. Keyword `mean` is added to the
inertia in these cases, e.g. with Euclidean and Manhattan distances.
Inertia is based on squared index, and keyword `squared` is added to the
name of distance, unless data were square root transformed (argument
`sqrt.dist=TRUE`). If an additive constant was used with argument `add`,
`Lingoes` or `Cailliez adjusted` is added to the the name of inertia,
and the value of the constant is printed.

## See also

[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md),
[`dist`](https://rdrr.io/r/stats/dist.html),
[`cmdscale`](https://rdrr.io/r/stats/cmdscale.html),
[`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md)
for underlying and related functions. Function
[`sppscores`](https://vegandevs.github.io/vegan/reference/sppscores.md)
can add species scores or replace existing species scores.

The function returns similar result object as
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) (see
[`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.md)).
This section for
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) gives a more
complete list of functions that can be used to access and analyse dbRDA
results.

## Examples

``` r
data(varespec, varechem)
## dbrda
dbrda(varespec ~ N + P + K + Condition(Al), varechem, dist="bray")
#> 
#> Call: dbrda(formula = varespec ~ N + P + K + Condition(Al), data =
#> varechem, distance = "bray")
#> 
#>               Inertia Proportion Rank RealDims
#> Total          4.5444     1.0000              
#> Conditional    0.9726     0.2140    1         
#> Constrained    0.9731     0.2141    3        3
#> Unconstrained  2.5987     0.5718   19       13
#> 
#> Inertia is squared Bray distance
#> 
#> Eigenvalues for constrained axes:
#> dbRDA1 dbRDA2 dbRDA3 
#> 0.5362 0.3198 0.1171 
#> 
#> Eigenvalues for unconstrained axes:
#>   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#> 0.9054 0.5070 0.3336 0.2581 0.2027 0.1605 0.1221 0.0825 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
## avoid negative eigenvalues with sqrt distances
dbrda(varespec ~ N + P + K + Condition(Al), varechem, dist="bray",
     sqrt.dist = TRUE)
#> 
#> Call: dbrda(formula = varespec ~ N + P + K + Condition(Al), data =
#> varechem, distance = "bray", sqrt.dist = TRUE)
#> 
#>               Inertia Proportion Rank
#> Total          6.9500     1.0000     
#> Conditional    0.9535     0.1372    1
#> Constrained    1.2267     0.1765    3
#> Unconstrained  4.7698     0.6863   19
#> 
#> Inertia is Bray distance
#> 
#> Eigenvalues for constrained axes:
#> dbRDA1 dbRDA2 dbRDA3 
#> 0.5817 0.4086 0.2365 
#> 
#> Eigenvalues for unconstrained axes:
#>   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#> 0.9680 0.6100 0.4469 0.3837 0.3371 0.3012 0.2558 0.2010 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
## avoid negative eigenvalues also with Jaccard distances
(m <- dbrda(varespec ~ N + P + K + Condition(Al), varechem, dist="jaccard"))
#> 
#> Call: dbrda(formula = varespec ~ N + P + K + Condition(Al), data =
#> varechem, distance = "jaccard")
#> 
#>               Inertia Proportion Rank
#> Total          6.5044     1.0000     
#> Conditional    1.0330     0.1588    1
#> Constrained    1.2068     0.1855    3
#> Unconstrained  4.2646     0.6557   19
#> 
#> Inertia is squared Jaccard distance
#> 
#> Eigenvalues for constrained axes:
#> dbRDA1 dbRDA2 dbRDA3 
#> 0.5992 0.3994 0.2082 
#> 
#> Eigenvalues for unconstrained axes:
#>   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#> 1.0388 0.6441 0.4518 0.3759 0.3239 0.2785 0.2279 0.1644 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
## add species scores
sppscores(m) <- wisconsin(varespec)
## pco
pco(varespec, dist = "bray", sqrt.dist = TRUE)
#> 
#> Call: pco(X = varespec, dist = "bray", sqrt.dist = TRUE)
#> 
#>               Inertia Rank
#> Total            6.95     
#> Unconstrained    6.95   23
#> 
#> Inertia is Bray distance
#> 
#> Eigenvalues for unconstrained axes:
#>   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#> 1.6348 1.1428 0.5658 0.4780 0.3737 0.3716 0.3074 0.2665 
#> (Showing 8 of 23 unconstrained eigenvalues)
#> 
```
