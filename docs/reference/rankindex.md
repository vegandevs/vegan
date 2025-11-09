# Compares Dissimilarity Indices for Gradient Detection

Rank correlations between dissimilarity indices and gradient separation.

## Usage

``` r
rankindex(grad, veg, indices = c("euc", "man", "gow", "bra", "kul"),
          stepacross = FALSE, method = "spearman", 
    metric = c("euclidean", "mahalanobis", "manhattan", "gower"),
    ...)
```

## Arguments

- grad:

  The gradient variable or matrix.

- veg:

  The community data matrix.

- indices:

  Dissimilarity indices compared, partial matches to alternatives in
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md).
  Alternatively, it can be a (named) list of functions returning objects
  of class 'dist'.

- stepacross:

  Use
  [`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md)
  to find a shorter path dissimilarity. The dissimilarities for site
  pairs with no shared species are set `NA` using
  [`no.shared`](https://vegandevs.github.io/vegan/reference/distconnected.md)
  so that indices with no fixed upper limit can also be analysed.

- method:

  Correlation method used.

- metric:

  Metric to evaluate the gradient separation. See Details.

- ...:

  Other parameters to
  [`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md).

## Details

A good dissimilarity index for multidimensional scaling should have a
high rank-order similarity with gradient separation. The function
compares most indices in
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md)
against gradient separation using rank correlation coefficients in
[`cor`](https://rdrr.io/r/stats/cor.html). The gradient separation
between each point is assessed using given `metric`. The default is to
use Euclidean distance of continuous variables scaled to unit variance,
or to use Gower metric for mixed data using function
[`daisy`](https://rdrr.io/pkg/cluster/man/daisy.html) when `grad` has
factors. The other alternatives are Mahalanabis distances which are
based on `grad` matrix scaled so that columns are orthogonal
(uncorrelated) and have unit variance, or Manhattan distances of `grad`
variables scaled to unit range.

The `indices` argument can accept any dissimilarity indices besides the
ones calculated by the
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md)
function. For this, the argument value should be a (possibly named) list
of functions. Each function must return a valid 'dist' object with
dissimilarities, similarities are not accepted and should be converted
into dissimilarities beforehand.

## Value

Returns a named vector of rank correlations.

## References

Faith, F.P., Minchin, P.R. and Belbin, L. (1987). Compositional
dissimilarity as a robust measure of ecological distance. *Vegetatio*
69, 57-68.

## Author

Jari Oksanen, with additions from Peter Solymos

## Note

There are several problems in using rank correlation coefficients.
Typically there are very many ties when \\n(n-1)/2\\ gradient separation
values are derived from just \\n\\ observations. Due to floating point
arithmetics, many tied values differ by machine epsilon and are
arbitrarily ranked differently by
[`rank`](https://rdrr.io/r/base/rank.html) used in
[`cor.test`](https://rdrr.io/r/stats/cor.test.html). Two indices which
are identical with certain transformation or standardization may differ
slightly (magnitude \\10^{-15}\\) and this may lead into third or fourth
decimal instability in rank correlations. Small differences in rank
correlations should not be taken too seriously. Probably this method
should be replaced with a sounder method, but I do not yet know which...
You may experiment with
[`mantel`](https://vegandevs.github.io/vegan/reference/mantel.md),
[`anosim`](https://vegandevs.github.io/vegan/reference/anosim.md) or
even
[`protest`](https://vegandevs.github.io/vegan/reference/procrustes.md).

Earlier version of this function used `method = "kendall"`, but that is
far too slow in large data sets.

The functions returning dissimilarity objects should be self contained,
because the `...` argument passes additional parameters to
[`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md)
and not to the functions supplied via the `indices` argument.

## See also

[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.md),
[`stepacross`](https://vegandevs.github.io/vegan/reference/stepacross.md),
[`no.shared`](https://vegandevs.github.io/vegan/reference/distconnected.md),
[`monoMDS`](https://vegandevs.github.io/vegan/reference/monoMDS.md),
[`cor`](https://rdrr.io/r/stats/cor.html),
[`Machine`](https://rdrr.io/r/base/base-defunct.html), and for
alternatives
[`anosim`](https://vegandevs.github.io/vegan/reference/anosim.md),
[`mantel`](https://vegandevs.github.io/vegan/reference/mantel.md) and
[`protest`](https://vegandevs.github.io/vegan/reference/procrustes.md).

## Examples

``` r
data(varespec)
data(varechem)
## The variables are automatically scaled
rankindex(varechem, varespec)
#>       euc       man       gow       bra       kul 
#> 0.2396330 0.2735087 0.2288358 0.2837910 0.2839834 
rankindex(varechem, wisconsin(varespec))
#>       euc       man       gow       bra       kul 
#> 0.4200990 0.4215642 0.3708606 0.4215642 0.4215642 
## Using non vegdist indices as functions
funs <- list(Manhattan=function(x) dist(x, "manhattan"),
    Gower=function(x) cluster:::daisy(x, "gower"),
    Ochiai=function(x) designdist(x, "1-J/sqrt(A*B)"))
rankindex(scale(varechem), varespec, funs)
#> Manhattan     Gower    Ochiai 
#> 0.2735087 0.2288358 0.1696862 
```
