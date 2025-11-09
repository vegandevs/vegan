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
  formula will add a unique indentifier to rows and constant for the
  rows to always produce estimates of row-level alpha and overall gamma
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
#>         statistic      SES    mean    2.5%     50%   97.5% Pr(sim.)  
#> alpha.1    11.235  -89.880 14.0911 14.0486 14.0873 14.1307     0.05 *
#> gamma      12.006 -289.268 14.1352 14.1234 14.1354 14.1459     0.05 *
#> beta.1      1.071   32.932  1.0031  1.0005  1.0032  1.0062     0.05 *
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
#>         statistic      SES     mean     2.5%      50%   97.5% Pr(sim.)  
#> alpha.1    8.0555  -65.965 12.18726 12.08926 12.18293 12.2823     0.05 *
#> alpha.2   11.2353  -70.080 14.09261 14.03438 14.09102 14.1628     0.05 *
#> alpha.3   12.0064 -395.831 14.13692 14.12571 14.13864 14.1432     0.05 *
#> gamma     14.1603    0.000 14.16027 14.16027 14.16027 14.1603     1.00  
#> beta.1     1.3568   31.669  1.16018  1.15130  1.16058  1.1701     0.05 *
#> beta.2     1.0710   24.436  1.00315  0.99857  1.00325  1.0073     0.05 *
#> beta.3     1.1794  465.895  1.00165  1.00121  1.00153  1.0024     0.05 *
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
#> alpha.1  8.055481  -53.689 12.202000 12.084996 12.193911 12.3400     0.05 *
#> alpha.2 11.235261 -111.437 14.102338 14.067943 14.092888 14.1472     0.05 *
#> alpha.3 12.006443 -440.415 14.136561 14.127995 14.136807 14.1443     0.05 *
#> gamma   14.160271    0.000 14.160271 14.160271 14.160271 14.1603     1.00  
#> beta.1   0.078594   20.640  0.068205  0.067538  0.068084  0.0691     0.05 *
#> beta.2   0.535514   37.487  0.501216  0.499667  0.501459  0.5025     0.05 *
#> beta.3   0.589695  518.507  0.500839  0.500565  0.500830  0.5011     0.05 *
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
#> alpha.1    8.0555  -56.207 12.2004 12.0598 12.2083 12.3284     0.05 *
#> alpha.2   11.2353  -96.350 14.0825 14.0365 14.0817 14.1365     0.05 *
#> alpha.3   12.0064 -314.260 14.1342 14.1226 14.1363 14.1427     0.05 *
#> gamma     14.1603    0.000 14.1603 14.1603 14.1603 14.1603     1.00  
#> beta.1     1.7578   85.055  1.1607  1.1486  1.1599  1.1742     0.05 *
#> beta.2     1.2603  120.827  1.0055  1.0017  1.0056  1.0088     0.05 *
#> beta.3     1.1794  369.880  1.0018  1.0012  1.0017  1.0027     0.05 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
