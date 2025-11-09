# Add or Drop Single Terms to a Constrained Ordination Model

Compute all single terms that can be added to or dropped from a
constrained ordination model.

## Usage

``` r
# S3 method for class 'cca'
add1(object, scope, test = c("none", "permutation"),
    permutations = how(nperm=199), ...)
# S3 method for class 'cca'
drop1(object, scope, test = c("none", "permutation"), 
    permutations = how(nperm=199), ...)
```

## Arguments

- object:

  A constrained ordination object from
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) or
  [`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md).

- scope:

  A formula giving the terms to be considered for adding or dropping;
  see [`add1`](https://rdrr.io/r/stats/add1.html) for details.

- test:

  Should a permutation test be added using
  [`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md).

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required, or a permutation matrix where each
  row gives the permuted indices.

- ...:

  Other arguments passed to
  [`add1.default`](https://rdrr.io/r/stats/add1.html),
  [`drop1.default`](https://rdrr.io/r/stats/add1.html), and
  [`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md).

## Details

With argument `test = "none"` the functions will only call
[`add1.default`](https://rdrr.io/r/stats/add1.html) or
[`drop1.default`](https://rdrr.io/r/stats/add1.html). With argument
`test = "permutation"` the functions will add test results from
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md).
Function `drop1.cca` will call
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md)
with argument `by = "margin"`. Function `add1.cca` will implement a test
for single term additions that is not directly available in
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md).

Functions are used implicitly in
[`step`](https://rdrr.io/r/stats/step.html),
[`ordiR2step`](https://vegandevs.github.io/vegan/reference/ordistep.md)
and
[`ordistep`](https://vegandevs.github.io/vegan/reference/ordistep.md).
The
[`deviance.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)
and
[`deviance.rda`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)
used in [`step`](https://rdrr.io/r/stats/step.html) have no firm basis,
and setting argument `test = "permutation"` may help in getting useful
insight into validity of model building. Function
[`ordistep`](https://vegandevs.github.io/vegan/reference/ordistep.md)
calls alternately `drop1.cca` and `add1.cca` with argument
`test = "permutation"` and selects variables by their permutation
\\P\\-values. Meticulous use of `add1.cca` and `drop1.cca` will allow
more judicious model building.

The default number of `permutations` is set to a low value, because
permutation tests can take a long time. It should be sufficient to give
a impression on the significances of the terms, but higher values of
`permutations` should be used if \\P\\ values really are important.

## Value

Returns a similar object as [`add1`](https://rdrr.io/r/stats/add1.html)
and [`drop1`](https://rdrr.io/r/stats/add1.html).

## Author

Jari Oksanen

## See also

[`add1`](https://rdrr.io/r/stats/add1.html),
[`drop1`](https://rdrr.io/r/stats/add1.html) and
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md)
for basic methods. You probably need these functions with
[`step`](https://rdrr.io/r/stats/step.html) and
[`ordistep`](https://vegandevs.github.io/vegan/reference/ordistep.md).
Functions
[`deviance.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)
and
[`extractAIC.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md)
are used to produce the other arguments than test results in the output.

## Examples

``` r
data(dune)
data(dune.env)
## Automatic model building based on AIC but with permutation tests
step(cca(dune ~  1, dune.env), reformulate(names(dune.env)), test="perm")
#> Start:  AIC=87.66
#> dune ~ 1
#> 
#>              Df    AIC      F Pr(>F)   
#> + Moisture    3 86.608 2.2536  0.005 **
#> + Management  3 86.935 2.1307  0.005 **
#> + A1          1 87.411 2.1400  0.045 * 
#> <none>          87.657                 
#> + Manure      4 88.832 1.5251  0.030 * 
#> + Use         2 89.134 1.1431  0.200   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step:  AIC=86.61
#> dune ~ Moisture
#> 
#>              Df    AIC      F Pr(>F)   
#> <none>          86.608                 
#> + Management  3 86.813 1.4565  0.070 . 
#> + A1          1 86.992 1.2624  0.205   
#> + Use         2 87.259 1.2760  0.165   
#> + Manure      4 87.342 1.3143  0.080 . 
#> - Moisture    3 87.657 2.2536  0.005 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
## see ?ordistep to do the same, but based on permutation P-values
if (FALSE) { # \dontrun{
ordistep(cca(dune ~  1, dune.env), reformulate(names(dune.env)))
} # }
## Manual model building
## -- define the maximal model for scope
mbig <- rda(dune ~  ., dune.env)
#> 
#> Some constraints or conditions were aliased because they were redundant. This
#> can happen if terms are constant or linearly dependent (collinear): ‘Manure^4’
## -- define an empty model to start with
m0 <- rda(dune ~ 1, dune.env)
## -- manual selection and updating
add1(m0, scope=formula(mbig), test="perm")
#>            Df    AIC      F Pr(>F)   
#> <none>        89.620                 
#> A1          1 89.591 1.9217  0.015 * 
#> Moisture    3 87.707 2.5883  0.005 **
#> Management  3 87.082 2.8400  0.005 **
#> Use         2 91.032 1.1741  0.225   
#> Manure      4 89.232 1.9539  0.005 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
m0 <- update(m0, . ~ . + Management)
add1(m0, scope=formula(mbig), test="perm")
#>          Df    AIC      F Pr(>F)  
#> <none>      87.082                
#> A1        1 87.424 1.2965  0.270  
#> Moisture  3 85.567 1.9764  0.015 *
#> Use       2 88.284 1.0510  0.415  
#> Manure    3 87.517 1.3902  0.130  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
m0 <- update(m0, . ~ . + Moisture)
## -- included variables still significant?
drop1(m0, test="perm")
#>            Df    AIC      F Pr(>F)   
#> <none>        85.567                 
#> Management  3 87.707 2.1769  0.010 **
#> Moisture    3 87.082 1.9764  0.015 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
add1(m0, scope=formula(mbig), test="perm")
#>        Df    AIC      F Pr(>F)
#> <none>    85.567              
#> A1      1 86.220 0.8359   0.62
#> Use     2 86.842 0.8027   0.78
#> Manure  3 85.762 1.1225   0.34
```
