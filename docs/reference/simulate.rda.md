# Simulate Responses with Gaussian Error or Permuted Residuals for Constrained Ordination

Function simulates a response data frame so that it adds Gaussian error
to the fitted responses of Redundancy Analysis
([`rda`](https://vegandevs.github.io/vegan/reference/cca.md)),
Constrained Correspondence Analysis
([`cca`](https://vegandevs.github.io/vegan/reference/cca.md)) or
distance-based RDA
([`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)).
The function is a special case of generic
[`simulate`](https://rdrr.io/r/stats/simulate.html), and works similarly
as `simulate.lm`.

## Usage

``` r
# S3 method for class 'rda'
simulate(object, nsim = 1, seed = NULL, indx = NULL,
    rank = "full", correlated = FALSE, ...)
```

## Arguments

- object:

  an object representing a fitted
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md) or
  [`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)
  model.

- nsim:

  number of response matrices to be simulated. Only one dissimilarity
  matrix is returned for
  [`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md),
  and larger `nsim` is an error.

- seed:

  an object specifying if and how the random number generator should be
  initialized (‘seeded’). See
  [`simulate`](https://rdrr.io/r/stats/simulate.html) for details.

- indx:

  Index of residuals added to the fitted values, such as produced by
  [`shuffleSet`](https://rdrr.io/pkg/permute/man/shuffleSet.html) or
  [`sample`](https://rdrr.io/r/base/sample.html). The index can have
  duplicate entries so that bootstrapping is allowed. If `nsim` \\\>1\\,
  the output should be compliant with
  [`shuffleSet`](https://rdrr.io/pkg/permute/man/shuffleSet.html) with
  one line for each simulation. If `nsim` is missing, the number of rows
  of `indx` is used to define the number of simulations, but if `nsim`
  is given, it should match number of rows in `indx`. If null,
  parametric simulation is used and Gaussian error is added to the
  fitted values.

- rank:

  The rank of the constrained component: passed to
  [`predict.rda`](https://vegandevs.github.io/vegan/reference/predict.cca.md)
  or
  [`predict.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md).

- correlated:

  Are species regarded as correlated in parametric simulation or when
  `indx` is not given? If `correlated = TRUE`, multivariate Gaussian
  random error is generated, and if `FALSE`, Gaussian random error is
  generated separately for each species. The argument has no effect in
  [`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)
  which has no information on species.

- ...:

  additional optional arguments (ignored).

## Details

The implementation follows `"lm"` method of
[`simulate`](https://rdrr.io/r/stats/simulate.html), and adds Gaussian
(Normal) error to the fitted values
([`fitted.rda`](https://vegandevs.github.io/vegan/reference/predict.cca.md))
using function [`rnorm`](https://rdrr.io/r/stats/Normal.html) if
`correlated = FALSE` or
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) if
`correlated = TRUE`. The standard deviations
([`rnorm`](https://rdrr.io/r/stats/Normal.html)) or covariance matrices
for species ([`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html)) are
estimated from the residuals after fitting the constraints.
Alternatively, the function can take a permutation index that is used to
add permuted residuals (unconstrained component) to the fitted values.
Raw data are used in
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md). Internal
Chi-square transformed data are used in
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) within the
function, but the returned matrix is similar to the original input data.
The simulation is performed on internal metric scaling data in
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md), but
the function returns the Euclidean distances calculated from the
simulated data. The simulation uses only the real components, and the
imaginary dimensions are ignored.

## Value

If `nsim = 1`, returns a matrix or dissimilarities (in
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)) with
similar additional arguments on random number seed as
[`simulate`](https://rdrr.io/r/stats/simulate.html). If `nsim > 1`,
returns a similar array as returned by
[`simulate.nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
with similar attributes.

## Author

Jari Oksanen

## See also

[`simulate`](https://rdrr.io/r/stats/simulate.html) for the generic case
and for [`lm`](https://rdrr.io/r/stats/lm.html) objects, and
[`simulate.nullmodel`](https://vegandevs.github.io/vegan/reference/nullmodel.md)
for community null model simulation. Functions
[`fitted.rda`](https://vegandevs.github.io/vegan/reference/predict.cca.md)
and
[`fitted.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md)
return fitted values without the error component. See
[`rnorm`](https://rdrr.io/r/stats/Normal.html) and
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) (MASS package)
for simulating Gaussian random error.

## Examples

``` r
data(dune)
data(dune.env)
mod <- rda(dune ~  Moisture + Management, dune.env)
## One simulation
update(mod, simulate(mod) ~  .)
#> 
#> Call: rda(formula = simulate(mod) ~ Moisture + Management, data = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total         80.1851     1.0000     
#> Constrained   55.0230     0.6862    6
#> Unconstrained 25.1621     0.3138   13
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6 
#> 21.135 19.419  7.051  3.688  2.279  1.451 
#> 
#> Eigenvalues for unconstrained axes:
#>   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10  PC11  PC12  PC13 
#> 4.723 3.700 3.165 3.071 2.709 1.992 1.761 1.214 1.108 0.745 0.492 0.272 0.210 
#> 
## An impression of confidence regions of site scores
plot(mod, display="sites")
for (i in 1:5) lines(procrustes(mod, update(mod, simulate(mod) ~ .)), col="blue")

## Simulate a set of null communities with permutation of residuals
simulate(mod, indx = shuffleSet(nrow(dune), 99))
#> An object of class “simulate.rda” 
#> ‘simulate index’ method (abundance, non-sequential)
#> 20 x 30 matrix
#> Number of permuted matrices = 99 
#> 
```
