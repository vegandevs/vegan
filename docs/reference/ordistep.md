# Choose a Model by Permutation Tests in Constrained Ordination

Automatic stepwise model building for constrained ordination methods
([`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md),
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)). The
function `ordistep` is modelled after
[`step`](https://rdrr.io/r/stats/step.html) and can do forward, backward
and stepwise model selection using permutation tests. Function
`ordiR2step` performs forward model choice solely on adjusted \\R^2\\
and \\P\\-value.

## Usage

``` r
ordistep(object, scope, direction = c("both", "backward", "forward"),
   Pin = 0.05, Pout = 0.1, permutations = how(nperm = 199), steps = 50,
   trace = TRUE, ...)
ordiR2step(object, scope, Pin = 0.05, R2scope = TRUE,
   permutations = how(nperm = 499), trace = TRUE, R2permutations = 1000, ...)
```

## Arguments

- object:

  In `ordistep`, an ordination object inheriting from
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md) or
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md).

- scope:

  Defines the range of models examined in the stepwise search. This can
  be a list containing components `upper` and `lower`, both formulae. If
  it is a single item, it is interpreted the target scope, depending on
  the `direction`. If `direction` is `"forward"`, a single item is
  interpreted as the `upper` scope and the formula of the input `object`
  as the `lower` scope. See [`step`](https://rdrr.io/r/stats/step.html)
  for details. In `ordiR2step`, this defines the upper scope; it can
  also be an ordination object from with the model is extracted.

- direction:

  The mode of stepwise search, can be one of `"both"`, `"backward"`, or
  `"forward"`, with a default of `"both"`. If the `scope` argument is
  missing, the default for `direction` is `"backward"` in `ordistep`
  (and `ordiR2step` does not have this argument, but only works
  forward).

- Pin, Pout:

  Limits of permutation \\P\\-values for adding (`Pin`) a term to the
  model, or dropping (`Pout`) from the model. Term is added if \\P \le\\
  `Pin`, and removed if \\P \>\\ `Pout`.

- R2scope:

  Use adjusted \\R^2\\ as the stopping criterion: only models with lower
  adjusted \\R^2\\ than scope are accepted.

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required, or a permutation matrix where each
  row gives the permuted indices. This is passed to
  [`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md):
  see there for details.

- steps:

  Maximum number of iteration steps of dropping and adding terms.

- trace:

  If positive, information is printed during the model building. Larger
  values may give more information.

- R2permutations:

  Number of permutations used in the estimation of adjusted \\R^2\\ for
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md) using
  [`RsquareAdj`](https://vegandevs.github.io/vegan/reference/RsquareAdj.md).

- ...:

  Any additional arguments to
  [`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md)
  and
  [`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md).

## Details

The basic functions for model choice in constrained ordination are
[`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md)
and
[`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md).
With these functions, ordination models can be chosen with standard R
function [`step`](https://rdrr.io/r/stats/step.html) which bases the
term choice on AIC. AIC-like statistics for ordination are provided by
functions
[`deviance.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)
and
[`extractAIC.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)
(with similar functions for
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md)). Actually,
constrained ordination methods do not have AIC, and therefore the
[`step`](https://rdrr.io/r/stats/step.html) may not be trusted. This
function provides an alternative using permutation \\P\\-values.

Function `ordistep` defines the model, `scope` of models considered, and
`direction` of the procedure similarly as
[`step`](https://rdrr.io/r/stats/step.html). The function alternates
with `drop` and `add` steps and stops when the model was not changed
during one step. The `-` and `+` signs in the summary table indicate
which stage is performed. It is often sensible to have `Pout` \\\>\\
`Pin` in stepwise models to avoid cyclic adds and drops of single terms.

Function `ordiR2step` builds model forward so that it maximizes adjusted
\\R^2\\ (function
[`RsquareAdj`](https://vegandevs.github.io/vegan/reference/RsquareAdj.md))
at every step, and stopping when the adjusted \\R^2\\ starts to
decrease, or the adjusted \\R^2\\ of the `scope` is exceeded, or the
selected permutation \\P\\-value is exceeded (Blanchet et al. 2008). The
second criterion is ignored with option `R2scope = FALSE`, and the third
criterion can be ignored setting `Pin = 1` (or higher). The function
cannot be used if adjusted \\R^2\\ cannot be calculated. If the number
of predictors is higher than the number of observations, adjusted
\\R^2\\ is also unavailable. Such models can be analysed with
`R2scope = FALSE`, but the variable selection will stop if models become
overfitted and adjusted \\R^2\\ cannot be calculated, and the adjusted
\\R^2\\ will be reported as zero. The \\R^2\\ of
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) is based on
simulations (see
[`RsquareAdj`](https://vegandevs.github.io/vegan/reference/RsquareAdj.md))
and different runs of `ordiR2step` can give different results.

Functions `ordistep` (based on \\P\\ values) and `ordiR2step` (based on
adjusted \\R^2\\ and hence on eigenvalues) can select variables in
different order.

## Value

Functions return the selected model with one additional component,
`anova`, which contains brief information of steps taken. You can
suppress voluminous output during model building by setting
`trace = FALSE`, and find the summary of model history in the `anova`
item.

## References

Blanchet, F. G., Legendre, P. & Borcard, D. (2008) Forward selection of
explanatory variables. *Ecology* 89, 2623–2632.

## Author

Jari Oksanen

## See also

The function handles constrained ordination methods
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) and
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md). The
underlying functions are
[`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md)
and
[`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md),
and the function is modelled after standard
[`step`](https://rdrr.io/r/stats/step.html) (which also can be used
directly but uses AIC for model choice, see
[`extractAIC.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)).
Function `ordiR2step` builds upon
[`RsquareAdj`](https://vegandevs.github.io/vegan/reference/RsquareAdj.md).

## Examples

``` r
## See add1.cca for another example

### Dune data
data(dune)
data(dune.env)
mod0 <- rda(dune ~ 1, dune.env)  # Model with intercept only
mod1 <- rda(dune ~ ., dune.env)  # Model with all explanatory variables
#> 
#> Some constraints or conditions were aliased because they were redundant. This
#> can happen if terms are constant or linearly dependent (collinear): ‘Manure^4’

## With scope present, the default direction is "both"
mod <- ordistep(mod0, scope = formula(mod1))
#> 
#> Start: dune ~ 1 
#> 
#>              Df    AIC      F Pr(>F)   
#> + Management  3 87.082 2.8400  0.005 **
#> + Moisture    3 87.707 2.5883  0.005 **
#> + Manure      4 89.232 1.9539  0.005 **
#> + A1          1 89.591 1.9217  0.060 . 
#> + Use         2 91.032 1.1741  0.295   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Management 
#> 
#>              Df   AIC    F Pr(>F)   
#> - Management  3 89.62 2.84  0.005 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>            Df    AIC      F Pr(>F)   
#> + Moisture  3 85.567 1.9764  0.010 **
#> + Manure    3 87.517 1.3902  0.140   
#> + A1        1 87.424 1.2965  0.220   
#> + Use       2 88.284 1.0510  0.345   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Management + Moisture 
#> 
#>              Df    AIC      F Pr(>F)   
#> - Moisture    3 87.082 1.9764  0.015 * 
#> - Management  3 87.707 2.1769  0.005 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>          Df    AIC      F Pr(>F)
#> + Manure  3 85.762 1.1225  0.360
#> + A1      1 86.220 0.8359  0.525
#> + Use     2 86.842 0.8027  0.740
#> 
mod
#> 
#> Call: rda(formula = dune ~ Management + Moisture, data = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total         84.1237     1.0000     
#> Constrained   46.4249     0.5519    6
#> Unconstrained 37.6988     0.4481   13
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6 
#> 21.588 14.075  4.123  3.163  2.369  1.107 
#> 
#> Eigenvalues for unconstrained axes:
#>   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10  PC11  PC12  PC13 
#> 8.241 7.138 5.355 4.409 3.143 2.770 1.878 1.741 0.952 0.909 0.627 0.311 0.227 
#> 
## summary table of steps
mod$anova
#>              Df    AIC      F Pr(>F)   
#> + Management  3 87.082 2.8400  0.005 **
#> + Moisture    3 85.567 1.9764  0.010 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Example of ordistep, forward
ordistep(mod0, scope = formula(mod1), direction="forward")
#> 
#> Start: dune ~ 1 
#> 
#>              Df    AIC      F Pr(>F)   
#> + Moisture    3 87.707 2.5883  0.005 **
#> + Management  3 87.082 2.8400  0.010 **
#> + Manure      4 89.232 1.9539  0.035 * 
#> + A1          1 89.591 1.9217  0.055 . 
#> + Use         2 91.032 1.1741  0.225   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Moisture 
#> 
#>              Df    AIC      F Pr(>F)   
#> + Management  3 85.567 2.1769  0.010 **
#> + Manure      4 86.060 1.8598  0.020 * 
#> + Use         2 88.003 1.4245  0.145   
#> + A1          1 88.886 0.6286  0.825   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Moisture + Management 
#> 
#>          Df    AIC      F Pr(>F)
#> + Manure  3 85.762 1.1225  0.360
#> + A1      1 86.220 0.8359  0.645
#> + Use     2 86.842 0.8027  0.725
#> 
#> 
#> Call: rda(formula = dune ~ Moisture + Management, data = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total         84.1237     1.0000     
#> Constrained   46.4249     0.5519    6
#> Unconstrained 37.6988     0.4481   13
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6 
#> 21.588 14.075  4.123  3.163  2.369  1.107 
#> 
#> Eigenvalues for unconstrained axes:
#>   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10  PC11  PC12  PC13 
#> 8.241 7.138 5.355 4.409 3.143 2.770 1.878 1.741 0.952 0.909 0.627 0.311 0.227 
#> 

## Example of ordiR2step (always forward)
## stops because R2 of 'mod1' exceeded
ordiR2step(mod0, mod1)
#> Step: R2.adj= 0 
#> Call: dune ~ 1 
#>  
#>                 R2.adjusted
#> <All variables>  0.32508817
#> + Management     0.22512409
#> + Moisture       0.20050225
#> + Manure         0.16723149
#> + A1             0.04626579
#> + Use            0.01799755
#> <model>          0.00000000
#> 
#>              Df    AIC    F Pr(>F)   
#> + Management  3 87.082 2.84  0.002 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: R2.adj= 0.2251241 
#> Call: dune ~ Management 
#>  
#>                 R2.adjusted
#> + Moisture        0.3450334
#> <All variables>   0.3250882
#> + Manure          0.2779515
#> + A1              0.2392216
#> + Use             0.2300349
#> <model>           0.2251241
#> 
#> 
#> Call: rda(formula = dune ~ Management, data = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total         84.1237     1.0000     
#> Constrained   29.2307     0.3475    3
#> Unconstrained 54.8930     0.6525   16
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3 
#> 14.865 10.690  3.675 
#> 
#> Eigenvalues for unconstrained axes:
#>    PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
#> 15.270  8.428  6.899  5.675  3.988  3.121  2.588  2.380  1.818  1.376  0.995 
#>   PC12   PC13   PC14   PC15   PC16 
#>  0.785  0.661  0.467  0.283  0.159 
#> 
```
