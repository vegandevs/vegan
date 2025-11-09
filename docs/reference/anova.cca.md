# Permutation Test for Constrained Correspondence Analysis, Redundancy Analysis and Constrained Analysis of Principal Coordinates

The function performs an ANOVA like permutation test for Constrained
Correspondence Analysis
([`cca`](https://vegandevs.github.io/vegan/reference/cca.md)),
Redundancy Analysis
([`rda`](https://vegandevs.github.io/vegan/reference/cca.md)) or
distance-based Redundancy Analysis (dbRDA,
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md)) to
assess the significance of constraints.

## Usage

``` r
# S3 method for class 'cca'
anova(object, ..., permutations = how(nperm=999),
     by = NULL, model = c("reduced", "direct", "full"),
     parallel = getOption("mc.cores"), strata = NULL,
     cutoff = 1, scope = NULL)
# S3 method for class 'cca'
permutest(x, permutations = how(nperm = 99),
     model = c("reduced", "direct", "full"), by = NULL, first = FALSE,
     strata = NULL, parallel = getOption("mc.cores"), ...)
```

## Arguments

- object:

  One or several result objects from
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) or
  [`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md). If
  there are several result objects, they are compared against each other
  in the order they were supplied. For a single object, a test specified
  in `by` or an overall test is given.

- x:

  A single ordination result object.

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required, or a permutation matrix where each
  row gives the permuted indices.

- by:

  Setting `by = "axis"` will assess significance for each constrained
  axis, and setting `by = "terms"` will assess significance for each
  term (sequentially from first to last), and setting `by = "margin"`
  will assess the marginal effects of the terms (each marginal term
  analysed in a model with all other variables), and `by = "onedf"` will
  assess sequentially one-degree-of-freedom contrasts of split factors.

- model:

  Permutation model: `model="direct"` permutes community data,
  `model="reduced"` permutes residuals of the community data after
  Conditions (partial model), `model = "full"` permutes residuals after
  Conditions and Constraints.

- parallel:

  Use parallel processing with the given number of cores.

- strata:

  An integer vector or factor specifying the strata for permutation. If
  supplied, observations are permuted only within the specified strata.
  It is an error to use this when `permutations` is a matrix, or a
  [`how`](https://rdrr.io/pkg/permute/man/how.html) defines `blocks`.
  This is a legacy argument that will be deprecated in the future: use
  `permutations = how(..., blocks)` instead.

- cutoff:

  Only effective with `by="axis"` where stops permutations after an axis
  equals or exceeds the `cutoff` \\p\\-value.

- scope:

  Only effective with `by="margin"` where it can be used to select the
  marginal terms for testing. The default is to test all marginal terms
  in [`drop.scope`](https://rdrr.io/r/stats/factor.scope.html).

- first:

  Analyse only significance of the first axis.

- ...:

  Parameters passed to other functions. `anova.cca` passes all arguments
  to `permutest.cca`. In `anova` with `by = "axis"` you can use argument
  `cutoff` (defaults `1`) which stops permutations after exceeding the
  given level.

## Details

Functions `anova.cca` and `permutest.cca` implement ANOVA like
permutation tests for the joint effect of constraints in
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) or
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md).
Function `anova` is intended as a more user-friendly alternative to
`permutest` (that is the real workhorse).

Function `anova` can analyse a sequence of constrained ordination
models. The analysis is based on the differences in residual deviance in
permutations of nested models.

The default test is for the sum of all constrained eigenvalues. Setting
`first = TRUE` will perform a test for the first constrained eigenvalue.
Argument `first` can be set either in `anova.cca` or in `permutest.cca`.
It is also possible to perform significance tests for each axis or for
each term (constraining variable) using argument `by` in `anova.cca`.
Setting `by = "axis"` will perform separate significance tests for each
constrained axis. All previous constrained axes will be used as
conditions (“partialled out”) and a test for the first constrained
eigenvalues is performed (Legendre et al. 2011). You can stop
permutation tests after exceeding a given significance level with
argument `cutoff` to speed up calculations in large models. Setting
`by = "terms"` will perform separate significance test for each term
(constraining variable). The terms are assessed sequentially from first
to last, and the order of the terms will influence their significances.
Setting `by = "onedf"` will perform a similar sequential test for
one-degree-of-freedom effects, where multi-level factors are split in
their contrasts. Setting `by = "margin"` will perform separate
significance test for each marginal term in a model with all other
terms. The marginal test also accepts a `scope` argument for the
[`drop.scope`](https://rdrr.io/r/stats/factor.scope.html) which can be a
character vector of term labels that are analysed, or a fitted model of
lower scope. The marginal effects are also known as “Type III” effects,
but the current function only evaluates marginal terms. It will, for
instance, ignore main effects that are included in interaction terms. In
calculating pseudo-\\F\\, all terms are compared to the same residual of
the full model.

Community data are permuted with choice `model="direct"`, and residuals
after partial CCA/ RDA/ dbRDA with choice `model="reduced"` (default).
If there is no partial CCA/ RDA/ dbRDA stage, `model="reduced"` simply
permutes the data and is equivalent to `model="direct"`. The test
statistic is “pseudo-\\F\\”, which is the ratio of constrained and
unconstrained total Inertia (Chi-squares, variances or something
similar), each divided by their respective degrees of freedom. If there
are no conditions (“partial” terms), the sum of all eigenvalues remains
constant, so that pseudo-\\F\\ and eigenvalues would give equal results.
In partial CCA/ RDA/ dbRDA, the effect of conditioning variables
(“covariables”) is removed before permutation, and the total Chi-square
is not fixed, and test based on pseudo-\\F\\ would differ from the test
based on plain eigenvalues.

## Value

The function `anova.cca` calls `permutest.cca` and fills an
[`anova`](https://rdrr.io/r/stats/anova.html) table. Additional
attributes are `Random.seed` (the random seeds used), `control` (the
permutation design, see [how](https://rdrr.io/pkg/permute/man/how.html))
and `F.perm` (the permuted test statistics).

## References

Legendre, P. and Legendre, L. (2012). *Numerical Ecology*. 3rd English
ed. Elsevier.

Legendre, P., Oksanen, J. and ter Braak, C.J.F. (2011). Testing the
significance of canonical axes in redundancy analysis. *Methods in
Ecology and Evolution* 2, 269–277.

## Author

Jari Oksanen

## See also

`anova.cca`,
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) to get
something to analyse. Function
[`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md)
calls `anova.cca` with `by = "margin"`, and
[`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md) an
analysis for single terms additions, which can be used in automatic or
semiautomatic model building (see
[`deviance.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)).

## Examples

``` r
data(dune, dune.env)
mod <- cca(dune ~ Moisture + Management, dune.env)
## overall test
anova(mod)
#> Permutation test for cca under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: cca(formula = dune ~ Moisture + Management, data = dune.env)
#>          Df ChiSquare      F Pr(>F)   
#> Model     6    1.0024 1.9515  0.003 **
#> Residual 13    1.1129                 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## tests for individual terms
anova(mod, by="term")
#> Permutation test for cca under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: cca(formula = dune ~ Moisture + Management, data = dune.env)
#>            Df ChiSquare      F Pr(>F)    
#> Moisture    3   0.62831 2.4465  0.001 ***
#> Management  3   0.37407 1.4565  0.049 *  
#> Residual   13   1.11289                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(mod, by="margin")
#> Permutation test for cca under reduced model
#> Marginal effects of terms
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: cca(formula = dune ~ Moisture + Management, data = dune.env)
#>            Df ChiSquare      F Pr(>F)  
#> Moisture    3   0.39854 1.5518  0.036 *
#> Management  3   0.37407 1.4565  0.053 .
#> Residual   13   1.11289                
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## sequential test for contrasts
anova(mod, by = "onedf")
#> Permutation test for cca under reduced model
#> Sequential test for contrasts
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: cca(formula = dune ~ Moisture + Management, data = dune.env)
#>              Df ChiSquare      F Pr(>F)    
#> Moisture.L    1   0.41081 4.7988  0.001 ***
#> Moisture.Q    1   0.11261 1.3154  0.167    
#> Moisture.C    1   0.10489 1.2253  0.226    
#> ManagementHF  1   0.08849 1.0337  0.359    
#> ManagementNM  1   0.20326 2.3744  0.011 *  
#> ManagementSF  1   0.08231 0.9615  0.455    
#> Residual     13   1.11289                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## test for adding all environmental variables
anova(mod, cca(dune ~ ., dune.env))
#> 
#> Some constraints or conditions were aliased because they were redundant. This
#> can happen if terms are constant or linearly dependent (collinear): ‘Manure^4’
#> Permutation tests for cca under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model 1: dune ~ Moisture + Management
#> Model 2: dune ~ A1 + Moisture + Management + Use + Manure
#>   ResDf ResChiSquare Df ChiSquare      F Pr(>F)
#> 1    13       1.1129                           
#> 2     7       0.6121  6   0.50079 0.9545  0.532
```
