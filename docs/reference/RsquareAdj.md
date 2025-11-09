# Adjusted R-square

The functions finds the adjusted R-square.

## Usage

``` r
# Default S3 method
RsquareAdj(x, n, m, ...)
# S3 method for class 'rda'
RsquareAdj(x, ...)
# S3 method for class 'cca'
RsquareAdj(x, permutations = 1000, ...)
```

## Arguments

- x:

  Unadjusted R-squared or an object from which the terms for evaluation
  or adjusted R-squared can be found.

- n, m:

  Number of observations and number of degrees of freedom in the fitted
  model.

- permutations:

  Number of permutations to use when computing the adjusted R-squared
  for a cca. The permutations can be calculated in parallel by
  specifying the number of cores which is passed to
  [`permutest`](https://vegandevs.github.io/vegan/reference/anova.cca.md)

- ...:

  Other arguments (ignored) except in the case of cca in which these
  arguments are passed to
  [`permutest`](https://vegandevs.github.io/vegan/reference/anova.cca.md).

## Details

The default method finds the adjusted \\R^2\\ from the unadjusted
\\R^2\\, number of observations, and number of degrees of freedom in the
fitted model. The specific methods find this information from the fitted
result object. There are specific methods for
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) (also used
for distance-based RDA),
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`lm`](https://rdrr.io/r/stats/lm.html) and
[`glm`](https://rdrr.io/r/stats/glm.html). Adjusted, or even unadjusted,
\\R^2\\ may not be available in some cases, and then the functions will
return `NA`. \\R^2\\ values are available only for
[`gaussian`](https://rdrr.io/r/stats/family.html) models in
[`glm`](https://rdrr.io/r/stats/glm.html).

The adjusted, \\R^2\\ of `cca` is computed using a permutation approach
developed by Peres-Neto et al. (2006). By default 1000 permutations are
used.

## Value

The functions return a list of items `r.squared` and `adj.r.squared`.

## References

Legendre, P., Oksanen, J. and ter Braak, C.J.F. (2011). Testing the
significance of canonical axes in redundancy analysis. *Methods in
Ecology and Evolution* 2, 269–277.

Peres-Neto, P., P. Legendre, S. Dray and D. Borcard. 2006. Variation
partitioning of species data matrices: estimation and comparison of
fractions. *Ecology* 87, 2614–2625.

## See also

[`varpart`](https://vegandevs.github.io/vegan/reference/varpart.md) uses
`RsquareAdj`.

## Examples

``` r
data(mite)
data(mite.env)
## rda
m <- rda(decostand(mite, "hell") ~  ., mite.env)
RsquareAdj(m)
#> $r.squared
#> [1] 0.5265047
#> 
#> $adj.r.squared
#> [1] 0.4367038
#> 
## cca
m <- cca(decostand(mite, "hell") ~  ., mite.env)
RsquareAdj(m)
#> $r.squared
#> [1] 0.4471676
#> 
#> $adj.r.squared
#> [1] 0.3447024
#> 
## default method
RsquareAdj(0.8, 20, 5)
#> [1] 0.7285714
```
