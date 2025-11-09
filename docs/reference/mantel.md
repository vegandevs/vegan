# Mantel and Partial Mantel Tests for Dissimilarity Matrices

Function `mantel` finds the Mantel statistic as a matrix correlation
between two dissimilarity matrices, and function `mantel.partial` finds
the partial Mantel statistic as the partial matrix correlation between
three dissimilarity matrices. The significance of the statistic is
evaluated by permuting rows and columns of the first dissimilarity
matrix. Test is one-sided and only tests that distances are positively
correlated.

## Usage

``` r
mantel(xdis, ydis, method="pearson", permutations=999, strata = NULL,
    na.rm = FALSE, parallel = getOption("mc.cores"))
mantel.partial(xdis, ydis, zdis, method = "pearson", permutations = 999, 
    strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))
# S3 method for class 'mantel'
summary(object, ...)
```

## Arguments

- xdis, ydis, zdis:

  Distance object of class `"dist"` or symmetric square matrices of
  distances. Only the lower triangle of square matrices is used. The
  first object `xdis` will be permuted in permutation tests.

- method:

  Correlation method, as accepted by
  [`cor`](https://rdrr.io/r/stats/cor.html): `"pearson"`, `"spearman"`
  or `"kendall"`.

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required, or a permutation matrix where each
  row gives the permuted indices.

- strata:

  An integer vector or factor specifying the strata for permutation. If
  supplied, observations are permuted only within the specified strata.

- na.rm:

  Remove missing values in calculation of Mantel correlation. Use this
  option with care: Permutation tests can be biased, in particular if
  two matrices had missing values in matching positions.

- parallel:

  Number of parallel processes or a predefined socket cluster. With
  `parallel = 1` uses ordinary, non-parallel processing. The parallel
  processing is done with parallel package.

- object:

  Result object.

- ...:

  Arguments passed to
  [`summary.permustats`](https://vegandevs.github.io/vegan/reference/permustats.md)
  These include `alternative` to select the sidedness of the test.

## Details

Mantel statistic is simply a correlation between entries of two
dissimilarity matrices (some use cross products, but these are linearly
related). However, the significance cannot be directly assessed, because
there are \\N(N-1)/2\\ entries for just \\N\\ observations. Mantel
developed asymptotic test, but here we use permutations of \\N\\ rows
and columns of dissimilarity matrix. Only the first matrix (`xdist`)
will be permuted, and the second is kept constant. See
[`permutations`](https://vegandevs.github.io/vegan/reference/permutations.md)
for additional details on permutation tests in Vegan.

Partial Mantel statistic uses partial correlation conditioned on the
third matrix. Only the first matrix is permuted so that the correlation
structure between second and first matrices is kept constant. Although
`mantel.partial` silently accepts other methods than `"pearson"`,
partial correlations will probably be wrong with other methods.

Borcard & Legendre (2012) warn against using partial Mantel test and
recommend instead Mantel correlogram
([`mantel.correlog`](https://vegandevs.github.io/vegan/reference/mantel.correlog.md)).

The function uses [`cor`](https://rdrr.io/r/stats/cor.html), which
should accept alternatives `pearson` for product moment correlations and
`spearman` or `kendall` for rank correlations.

## Value

The function returns a list of class `mantel` with following components:

- Call :

  Function call.

- method :

  Correlation method used, as returned by
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html).

- statistic:

  The Mantel statistic.

- signif:

  Empirical significance level from permutations.

- perm:

  A vector of permuted values. The distribution of permuted values can
  be inspected with
  [`permustats`](https://vegandevs.github.io/vegan/reference/permustats.md)
  function.

- permutations:

  Number of permutations.

- control:

  A list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html).

## References

Borcard, D. & Legendre, P. (2012) Is the Mantel correlogram powerful
enough to be useful in ecological analysis? A simulation study.
*Ecology* 93: 1473-1481.

Legendre, P. and Legendre, L. (2012) *Numerical Ecology*. 3rd English
Edition. Elsevier.

## Author

Jari Oksanen

## See also

[`mantel.correlog`](https://vegandevs.github.io/vegan/reference/mantel.correlog.md).

## Examples

``` r
## Is vegetation related to environment?
data(varespec)
data(varechem)
veg.dist <- vegdist(varespec) # Bray-Curtis
env.dist <- vegdist(scale(varechem), "euclid")
mantel(veg.dist, env.dist)
#> 
#> Mantel statistic based on Pearson's product-moment correlation 
#> 
#> Call:
#> mantel(xdis = veg.dist, ydis = env.dist) 
#> 
#> Mantel statistic r: 0.3047 
#>       Significance: 0.001 
#> 
#> Upper quantiles of permutations (null model):
#>   90%   95% 97.5%   99% 
#> 0.112 0.142 0.169 0.201 
#> Permutation: free
#> Number of permutations: 999
#> 
mantel(veg.dist, env.dist, method="spear")
#> 
#> Mantel statistic based on Spearman's rank correlation rho 
#> 
#> Call:
#> mantel(xdis = veg.dist, ydis = env.dist, method = "spear") 
#> 
#> Mantel statistic r: 0.2838 
#>       Significance: 0.001 
#> 
#> Upper quantiles of permutations (null model):
#>   90%   95% 97.5%   99% 
#> 0.111 0.143 0.175 0.214 
#> Permutation: free
#> Number of permutations: 999
#> 
```
