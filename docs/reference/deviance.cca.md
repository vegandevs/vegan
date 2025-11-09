# Statistics Resembling Deviance and AIC for Constrained Ordination

The functions extract statistics that resemble deviance and AIC from the
result of constrained correspondence analysis
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) or
redundancy analysis
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md). These
functions are rarely needed directly, but they are called by
[`step`](https://rdrr.io/r/stats/step.html) in automatic model building.
Actually, [`cca`](https://vegandevs.github.io/vegan/reference/cca.md)
and [`rda`](https://vegandevs.github.io/vegan/reference/cca.md) do not
have [`AIC`](https://rdrr.io/r/stats/AIC.html) and these functions are
certainly wrong.

## Usage

``` r
# S3 method for class 'cca'
deviance(object, ...)

# S3 method for class 'cca'
extractAIC(fit, scale = 0, k = 2, ...)
```

## Arguments

- object:

  the result of a constrained ordination
  ([`cca`](https://vegandevs.github.io/vegan/reference/cca.md) or
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md)).

- fit:

  fitted model from constrained ordination.

- scale:

  optional numeric specifying the scale parameter of the model, see
  `scale` in [`step`](https://rdrr.io/r/stats/step.html).

- k:

  numeric specifying the "weight" of the *equivalent degrees of freedom*
  (=:`edf`) part in the AIC formula.

- ...:

  further arguments.

## Details

The functions find statistics that resemble
[`deviance`](https://rdrr.io/r/stats/deviance.html) and
[`AIC`](https://rdrr.io/r/stats/AIC.html) in constrained ordination.
Actually, constrained ordination methods do not have a log-Likelihood,
which means that they cannot have AIC and deviance. Therefore you should
not use these functions, and if you use them, you should not trust them.
If you use these functions, it remains as your responsibility to check
the adequacy of the result.

The deviance of
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) is equal to
the Chi-square of the residual data matrix after fitting the
constraints. The deviance of
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) is defined
as the residual sum of squares. The deviance function of `rda` is also
used for distance-based RDA
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md).
Function `extractAIC` mimics `extractAIC.lm` in translating deviance to
AIC.

There is little need to call these functions directly. However, they are
called implicitly in [`step`](https://rdrr.io/r/stats/step.html)
function used in automatic selection of constraining variables. You
should check the resulting model with some other criteria, because the
statistics used here are unfounded. In particular, the penalty `k` is
not properly defined, and the default `k = 2` is not justified
theoretically. If you have only continuous covariates, the `step`
function will base the model building on magnitude of eigenvalues, and
the value of `k` only influences the stopping point (but the variables
with the highest eigenvalues are not necessarily the most significant in
permutation tests in
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md)).
If you also have multi-class factors, the value of `k` will have a
capricious effect in model building. The
[`step`](https://rdrr.io/r/stats/step.html) function will pass arguments
to [`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md)
and
[`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md),
and setting `test = "permutation"` will provide permutation tests of
each deletion and addition which can help in judging the validity of the
model building.

## Value

The `deviance` functions return “deviance”, and `extractAIC` returns
effective degrees of freedom and “AIC”.

## References

Godínez-Domínguez, E. & Freire, J. (2003) Information-theoretic approach
for selection of spatial and temporal models of community organization.
*Marine Ecology Progress Series* **253**, 17–24.

## Author

Jari Oksanen

## Note

These functions are unfounded and untested and they should not be used
directly or implicitly. Moreover, usual caveats in using
[`step`](https://rdrr.io/r/stats/step.html) are very valid.

## See also

[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`step`](https://rdrr.io/r/stats/step.html),
[`extractAIC`](https://rdrr.io/r/stats/extractAIC.html),
[`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md),
[`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md).

## Examples

``` r
# The deviance of correspondence analysis equals Chi-square
data(dune)
data(dune.env)
chisq.test(dune)
#> Warning: Chi-squared approximation may be incorrect
#> 
#>  Pearson's Chi-squared test
#> 
#> data:  dune
#> X-squared = 1449, df = 551, p-value < 2.2e-16
#> 
deviance(cca(dune))
#> [1] 1448.956
# Stepwise selection (forward from an empty model "dune ~ 1")
ord <- cca(dune ~ ., dune.env)
#> 
#> Some constraints or conditions were aliased because they were redundant. This
#> can happen if terms are constant or linearly dependent (collinear): ‘Manure^4’
step(cca(dune ~ 1, dune.env), scope = formula(ord))
#> Start:  AIC=87.66
#> dune ~ 1
#> 
#>              Df    AIC
#> + Moisture    3 86.608
#> + Management  3 86.935
#> + A1          1 87.411
#> <none>          87.657
#> + Manure      4 88.832
#> + Use         2 89.134
#> 
#> Step:  AIC=86.61
#> dune ~ Moisture
#> 
#>              Df    AIC
#> <none>          86.608
#> + Management  3 86.813
#> + A1          1 86.992
#> + Use         2 87.259
#> + Manure      4 87.342
#> - Moisture    3 87.657
#> 
#> Call: cca(formula = dune ~ Moisture, data = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total          2.1153     1.0000     
#> Constrained    0.6283     0.2970    3
#> Unconstrained  1.4870     0.7030   16
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>   CCA1   CCA2   CCA3 
#> 0.4187 0.1330 0.0766 
#> 
#> Eigenvalues for unconstrained axes:
#>    CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8    CA9   CA10   CA11 
#> 0.4098 0.2259 0.1761 0.1234 0.1082 0.0908 0.0859 0.0609 0.0566 0.0467 0.0419 
#>   CA12   CA13   CA14   CA15   CA16 
#> 0.0201 0.0143 0.0099 0.0085 0.0080 
#> 
```
