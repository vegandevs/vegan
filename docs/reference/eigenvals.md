# Extract Eigenvalues from an Ordination Object

Function extracts eigenvalues from an object that has them. Many
multivariate methods return such objects.

## Usage

``` r
eigenvals(x, ...)
# S3 method for class 'cca'
eigenvals(x, model = c("all", "unconstrained", "constrained"),
          constrained = NULL, ...)
# S3 method for class 'decorana'
eigenvals(x, kind = c("additive", "axiswise", "decorana"),
           ...)
# S3 method for class 'eigenvals'
summary(object, ...)
```

## Arguments

- x:

  An object from which to extract eigenvalues.

- object:

  An `eigenvals` result object.

- model:

  Which eigenvalues to return for objects that inherit from class
  `"cca"` only.

- constrained:

  Return only constrained eigenvalues. Deprecated as of vegan 2.5-0. Use
  `model` instead.

- kind:

  Kind of eigenvalues returned for
  [`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md).
  Only `"additive"` eigenvalues can be used for reporting importances of
  components in `summary`. `"axiswise"` gives the non-additive
  eigenvalues, and `"decorana"` the decorana values (see
  [`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)
  for details).

- ...:

  Other arguments to the functions (usually ignored)

## Details

This is a generic function that has methods for
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md),
[`pcnm`](https://vegandevs.github.io/vegan/reference/pcnm.md),
[`prcomp`](https://rdrr.io/r/stats/prcomp.html),
[`princomp`](https://rdrr.io/r/stats/princomp.html), `dudi` (of ade4),
and `pca` and `pco` (of labdsv) result objects. The default method also
extracts eigenvalues if the result looks like being from
[`eigen`](https://rdrr.io/r/base/eigen.html) or
[`svd`](https://rdrr.io/r/base/svd.html). Functions
[`prcomp`](https://rdrr.io/r/stats/prcomp.html) and
[`princomp`](https://rdrr.io/r/stats/princomp.html) contain square roots
of eigenvalues that all called standard deviations, but `eigenvals`
function returns their squares. Function
[`svd`](https://rdrr.io/r/base/svd.html) contains singular values, but
function `eigenvals` returns their squares. For constrained ordination
methods [`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) and
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md) the
function returns the both constrained and unconstrained eigenvalues
concatenated in one vector, but the partial component will be ignored.
However, with argument `constrained = TRUE` only constrained eigenvalues
are returned.

The `summary` of `eigenvals` result returns eigenvalues, proportion
explained and cumulative proportion explained. The result object can
have some negative eigenvalues
([`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md),
[`pcnm`](https://vegandevs.github.io/vegan/reference/pcnm.md)) which
correspond to imaginary axes of Euclidean mapping of non-Euclidean
distances (Gower 1985). In these case real axes (corresponding to
positive eigenvalues) will "explain" proportion \>1 of total variation,
and negative eigenvalues bring the cumulative proportion to 1.
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md) will
only find the positive eigenvalues and only these are used in finding
proportions. For
[`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)
the importances and cumulative proportions are only reported for
`kind = "additive"`, because other alternatives do not add up to total
inertia of the input data.

## Value

An object of class `"eigenvals"`, which is a vector of eigenvalues.

The `summary` method returns an object of class `"summary.eigenvals"`,
which is a matrix.

## Author

Jari Oksanen.

## References

Gower, J. C. (1985). Properties of Euclidean and non-Euclidean distance
matrices. *Linear Algebra and its Applications* 67, 81–97.

## See also

[`eigen`](https://rdrr.io/r/base/eigen.html),
[`svd`](https://rdrr.io/r/base/svd.html),
[`prcomp`](https://rdrr.io/r/stats/prcomp.html),
[`princomp`](https://rdrr.io/r/stats/princomp.html),
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md),
[`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md),
[`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.md).

## Examples

``` r
data(varespec)
data(varechem)
mod <- cca(varespec ~ Al + P + K, varechem)
ev <- eigenvals(mod)
ev
#>      CCA1      CCA2      CCA3       CA1       CA2       CA3       CA4       CA5 
#> 0.3615566 0.1699600 0.1126167 0.3500372 0.2200788 0.1850741 0.1551179 0.1351054 
#>       CA6       CA7       CA8       CA9      CA10      CA11      CA12      CA13 
#> 0.1002670 0.0772991 0.0536938 0.0365603 0.0350887 0.0282291 0.0170651 0.0122474 
#>      CA14      CA15      CA16      CA17      CA18      CA19      CA20 
#> 0.0101910 0.0094701 0.0055090 0.0030529 0.0025118 0.0019485 0.0005178 
summary(ev)
#> Importance of components:
#>                         CCA1    CCA2    CCA3    CA1    CA2     CA3     CA4
#> Eigenvalue            0.3616 0.16996 0.11262 0.3500 0.2201 0.18507 0.15512
#> Proportion Explained  0.1736 0.08159 0.05406 0.1680 0.1056 0.08884 0.07446
#> Cumulative Proportion 0.1736 0.25514 0.30920 0.4772 0.5829 0.67172 0.74618
#>                           CA5     CA6     CA7     CA8     CA9    CA10    CA11
#> Eigenvalue            0.13511 0.10027 0.07730 0.05369 0.03656 0.03509 0.02823
#> Proportion Explained  0.06485 0.04813 0.03711 0.02577 0.01755 0.01684 0.01355
#> Cumulative Proportion 0.81104 0.85917 0.89627 0.92205 0.93960 0.95644 0.96999
#>                           CA12     CA13     CA14     CA15     CA16     CA17
#> Eigenvalue            0.017065 0.012247 0.010191 0.009470 0.005509 0.003053
#> Proportion Explained  0.008192 0.005879 0.004892 0.004546 0.002644 0.001465
#> Cumulative Proportion 0.978183 0.984062 0.988954 0.993500 0.996145 0.997610
#>                           CA18      CA19      CA20
#> Eigenvalue            0.002512 0.0019485 0.0005178
#> Proportion Explained  0.001206 0.0009353 0.0002486
#> Cumulative Proportion 0.998816 0.9997514 1.0000000

## choose which eignevalues to return
eigenvals(mod, model = "unconstrained")
#>       CA1       CA2       CA3       CA4       CA5       CA6       CA7       CA8 
#> 0.3500372 0.2200788 0.1850741 0.1551179 0.1351054 0.1002670 0.0772991 0.0536938 
#>       CA9      CA10      CA11      CA12      CA13      CA14      CA15      CA16 
#> 0.0365603 0.0350887 0.0282291 0.0170651 0.0122474 0.0101910 0.0094701 0.0055090 
#>      CA17      CA18      CA19      CA20 
#> 0.0030529 0.0025118 0.0019485 0.0005178 
```
