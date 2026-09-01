# Multiplicative Diversity Partitioning

In multiplicative diversity partitioning, mean values of alpha diversity
at lower levels of a sampling hierarchy are compared to the total
diversity in the entire data set or the pooled samples (gamma
diversity).

## Usage

``` r
multipart(...)
# Default S3 method
multipart(y, x, index=c("renyi", "tsallis"), scales = 1,
    global = FALSE, relative = FALSE, nsimul=99, method = "r2dtable", ...)
# S3 method for class 'formula'
multipart(formula, data, index=c("renyi", "tsallis"), scales = 1,
    global = FALSE, relative = FALSE, nsimul=99, method = "r2dtable", ...)
```

## Arguments

- y:

  A community matrix.

- x:

  A matrix with same number of rows as in `y`, columns coding the levels
  of sampling hierarchy. The number of groups within the hierarchy must
  decrease from left to right. If `x` is missing, two levels are
  assumed: each row is a group in the first level, and all rows are in
  the same group in the second level.

- formula:

  A two sided model formula in the form `y ~ x`, where `y` is the
  community data matrix with samples as rows and species as column.
  Right hand side (`x`) must be grouping variable(s) referring to levels
  of sampling hierarchy, terms from right to left will be treated as
  nested (first column is the lowest, last is the highest level). The
  formula will add a unique identifier to rows and constant for the rows
  to always produce estimates of row-level alpha and overall gamma
  diversities. You must use non-formula interface to avoid this
  behaviour. Interaction terms are not allowed.

- data:

  A data frame where to look for variables defined in the right hand
  side of `formula`. If missing, variables are looked in the global
  environment.

- index:

  Character, the entropy index to be calculated (see Details).

- relative:

  Logical, if `TRUE` then beta diversity is standardized by its maximum
  (see Details).

- scales:

  Numeric, of length 1, the order of the generalized diversity index to
  be used.

- global:

  Logical, indicates the calculation of beta diversity values, see
  Details.

- nsimul:

  Number of permutations to use. If `nsimul = 0`, only the `FUN`
  argument is evaluated. It is thus possible to reuse the statistic
  values without a null model.

- method:

  Null model method: either a name (character string) of a method
  defined in
  [`make.commsim`](https://vegandevs.github.io/vegan/reference/commsim.md)
  or a
  [`commsim`](https://vegandevs.github.io/vegan/reference/commsim.md)
  function. The default `"r2dtable"` keeps row sums and column sums
  fixed. See
  [`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md)
  for Details and Examples.

- ...:

  Other arguments passed to
  [`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md),
  i.e. `method`, `thin` or `burnin`.

## Details

Multiplicative diversity partitioning is based on Whittaker's (1972)
ideas, that has recently been generalised to one parametric diversity
families (i.e. Rényi and Tsallis) by Jost (2006, 2007). Jost recommends
to use the numbers equivalents (Hill numbers), instead of pure
diversities, and proofs, that this satisfies the multiplicative
partitioning requirements.

The current implementation of `multipart` calculates Hill numbers based
on the functions
[`renyi`](https://vegandevs.github.io/vegan/reference/renyi.md) and
[`tsallis`](https://vegandevs.github.io/vegan/reference/tsallis.md)
(provided as `index` argument). If values for more than one `scales` are
desired, it should be done in separate runs, because it adds extra
dimensionality to the implementation, which has not been resolved
efficiently.

Alpha diversities are then the averages of these Hill numbers for each
hierarchy levels, the global gamma diversity is the alpha value
calculated for the highest hierarchy level. When `global = TRUE`, beta
is calculated relative to the global gamma value: \$\$\beta_i = \gamma /
\alpha\_{i}\$\$ when `global = FALSE`, beta is calculated relative to
local gamma values (local gamma means the diversity calculated for a
particular cluster based on the pooled abundance vector): \$\$\beta_ij =
\alpha\_{(i+1)j} / mean(\alpha\_{ij})\$\$ where \\j\\ is a particular
cluster at hierarchy level \\i\\. Then beta diversity value for level
\\i\\ is the mean of the beta values of the clusters at that level,
\\\beta\_{i} = mean(\beta\_{ij})\\.

If `relative = TRUE`, the respective beta diversity values are
standardized by their maximum possible values (\\mean(\beta\_{ij}) /
\beta\_{max,ij}\\) given as \\\beta\_{max,ij} = n\_{j}\\ (the number of
lower level units in a given cluster \\j\\).

The expected diversity components are calculated `nsimul` times by
individual based randomization of the community data matrix. This is
done by the `"r2dtable"` method in
[`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md) by
default.

## Value

An object of class `"multipart"` with same structure as `"oecosimu"`
objects.

## References

Jost, L. (2006). Entropy and diversity. *Oikos*, **113**, 363–375.

Jost, L. (2007). Partitioning diversity into independent alpha and beta
components. *Ecology*, **88**, 2427–2439.

Whittaker, R. (1972). Evolution and measurement of species diversity.
*Taxon*, **21**, 213–251.

## Author

Péter Sólymos, <solymos@ualberta.ca>

## See also

See [`adipart`](https://vegandevs.github.io/vegan/reference/adipart.md)
for additive diversity partitioning,
[`hiersimu`](https://vegandevs.github.io/vegan/reference/adipart.md) for
hierarchical null model testing and
[`oecosimu`](https://vegandevs.github.io/vegan/reference/oecosimu.md)
for permutation settings and calculating \\p\\-values.

## Examples

``` r
## NOTE: 'nsimul' argument usually needs to be >= 99
## here much lower value is used for demonstration

data(mite)
data(mite.xy)
data(mite.env)
## Function to get equal area partitions of the mite data
cutter <- function (x, cut = seq(0, 10, by = 2.5)) {
    out <- rep(1, length(x))
    for (i in 2:(length(cut) - 1))
        out[which(x > cut[i] & x <= cut[(i + 1)])] <- i
    return(out)}
## The hierarchy of sample aggregation
levsm <- with(mite.xy, data.frame(
    l2=cutter(y, cut = seq(0, 10, by = 2.5)),
    l3=cutter(y, cut = seq(0, 10, by = 5))))
## Multiplicative diversity partitioning
multipart(mite, levsm, index="renyi", scales=1, nsimul=19)
#> multipart object
#> 
#> Call: multipart(y = mite, x = levsm, index = "renyi", scales = 1,
#> nsimul = 19)
#> 
#> nullmodel method ‘r2dtable’ with 19 simulations
#> options:  index renyi, scales 1, global FALSE
#> alternative hypothesis: statistic is less or greater than simulated values
#> 
#>         statistic      SES     mean     2.5%      50%   97.5% Pr(sim.)  
#> alpha.1    11.235  -82.989 14.08849 14.05444 14.07625 14.1661     0.05 *
#> gamma      12.006 -354.946 14.13614 14.12885 14.13494 14.1489     0.05 *
#> beta.1      1.071   28.451  1.00339  0.99804  1.00422  1.0058     0.05 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
multipart(mite ~ l2 + l3, levsm, index="renyi", scales=1, nsimul=19)
#> multipart object
#> 
#> Call: multipart(formula = mite ~ l2 + l3, data = levsm, index =
#> "renyi", scales = 1, nsimul = 19)
#> 
#> nullmodel method ‘r2dtable’ with 19 simulations
#> options:  index renyi, scales 1, global FALSE
#> alternative hypothesis: statistic is less or greater than simulated values
#> 
#>         statistic      SES    mean    2.5%     50%   97.5% Pr(sim.)  
#> alpha.1    8.0555  -59.908 12.1842 12.0956 12.1664 12.3348     0.05 *
#> alpha.2   11.2353  -82.190 14.0806 14.0265 14.0731 14.1400     0.05 *
#> alpha.3   12.0064 -324.898 14.1369 14.1249 14.1379 14.1464     0.05 *
#> gamma     14.1603    0.000 14.1603 14.1603 14.1603 14.1603     1.00  
#> beta.1     1.3568   26.899  1.1588  1.1446  1.1598  1.1705     0.05 *
#> beta.2     1.0710   28.310  1.0040  1.0002  1.0038  1.0078     0.05 *
#> beta.3     1.1794  382.497  1.0017  1.0010  1.0016  1.0025     0.05 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
multipart(mite ~ ., levsm, index="renyi", scales=1, nsimul=19, relative=TRUE)
#> multipart object
#> 
#> Call: multipart(formula = mite ~ ., data = levsm, index = "renyi",
#> scales = 1, relative = TRUE, nsimul = 19)
#> 
#> nullmodel method ‘r2dtable’ with 19 simulations
#> options:  index renyi, scales 1, global FALSE
#> alternative hypothesis: statistic is less or greater than simulated values
#> 
#>         statistic      SES      mean      2.5%       50%   97.5% Pr(sim.)  
#> alpha.1  8.055481  -58.630 12.188486 12.073974 12.180222 12.3118     0.05 *
#> alpha.2 11.235261  -88.711 14.089992 14.032972 14.085038 14.1477     0.05 *
#> alpha.3 12.006443 -337.258 14.136819 14.121477 14.138647 14.1439     0.05 *
#> gamma   14.160271    0.000 14.160271 14.160271 14.160271 14.1603     1.00  
#> beta.1   0.078594   22.888  0.068255  0.067282  0.068378  0.0688     0.05 *
#> beta.2   0.535514   30.034  0.501670  0.499752  0.501848  0.5037     0.05 *
#> beta.3   0.589695  396.846  0.500830  0.500578  0.500765  0.5014     0.05 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
multipart(mite ~ ., levsm, index="renyi", scales=1, nsimul=19, global=TRUE)
#> multipart object
#> 
#> Call: multipart(formula = mite ~ ., data = levsm, index = "renyi",
#> scales = 1, global = TRUE, nsimul = 19)
#> 
#> nullmodel method ‘r2dtable’ with 19 simulations
#> options:  index renyi, scales 1, global TRUE
#> alternative hypothesis: statistic is less or greater than simulated values
#> 
#>         statistic      SES    mean    2.5%     50%   97.5% Pr(sim.)  
#> alpha.1    8.0555  -64.055 12.1757 12.0786 12.1665 12.2886     0.05 *
#> alpha.2   11.2353  -94.848 14.0886 14.0355 14.0903 14.1396     0.05 *
#> alpha.3   12.0064 -349.235 14.1354 14.1242 14.1368 14.1455     0.05 *
#> gamma     14.1603    0.000 14.1603 14.1603 14.1603 14.1603     1.00  
#> beta.1     1.7578   96.885  1.1630  1.1523  1.1639  1.1723     0.05 *
#> beta.2     1.2603  118.875  1.0051  1.0015  1.0050  1.0089     0.05 *
#> beta.3     1.1794  411.140  1.0018  1.0010  1.0017  1.0026     0.05 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
