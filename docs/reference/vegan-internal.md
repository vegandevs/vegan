# Internal vegan functions

Internal vegan functions that are not intended to be called directly,
but only within other functions.

## Usage

``` r
ordiParseFormula(formula, data, xlev = NULL,  na.action = na.fail,
    subset = NULL, X)
ordiTerminfo(d, data)
ordiNAexclude(x, excluded)
ordiNApredict(omit, x)
ordiArgAbsorber(..., shrink, origin, scaling, triangular,
                display, choices, const, truemean, optimize, arrows, FUN)
centroids.cca(x, mf, wt)
getPermuteMatrix(perm, N, strata = NULL)
howHead(x, ...)
pasteCall(call, prefix = "Call:")
veganCovEllipse(cov, center = c(0, 0), scale = 1, npoints = 100)
veganMahatrans(x, s2, tol = sqrt(.Machine$double.eps), na.rm = FALSE)
hierParseFormula(formula, data)
GowerDblcen(x, na.rm = TRUE)
addLingoes(d)
addCailliez(d)
```

## Details

The description here is only intended for vegan developers: these
functions are not intended for users, but they only should be used
within functions. In general, these functions are not exported to the
namespace, but you must use [`get`](https://rdrr.io/r/base/get.html) or
[`:::`](https://rdrr.io/r/base/ns-dblcolon.html) to directly call these
functions.

`ordiParseFormula` returns a list of three matrices (dependent
variables, and
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) of
constraints and conditions, possibly `NULL`) needed in constrained
ordination. Argument `xlev` is passed to
[`model.frame`](https://rdrr.io/r/stats/model.frame.html). If the
left-hand-side was already evaluated in calling code, it can be given as
argument `X` and will not be re-evaluated. `ordiTermInfo` finds the term
information for constrained ordination as described in
[`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.md).
`ordiNAexclude` implements `na.action = na.exclude` for constrained
ordination finding WA scores of CCA components and site scores of
unconstrained component from `excluded` rows of observations. Function
`ordiNApredict` pads the result object with these or with WA scores
similarly as [`napredict`](https://rdrr.io/r/stats/nafns.html).

`ordiArgAbsorber` absorbs arguments of
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md)
function of vegan so that these do not cause superfluous warnings in
graphical function `FUN`. If you implement `scores` functions with new
arguments, you should update `ordiArgAbsorber`.

`centroids.cca` finds the weighted centroids of variables.

`getPermuteMatrix` interprets user input and returns a permutation
matrix where each row gives indices of observations for a permutation.
The input `perm` can be a single number for the number of simple
permutations, a result of
[`how`](https://rdrr.io/pkg/permute/man/how.html) defining a permutation
scheme or a permutation matrix.

`howHead` formats the permutation scheme of
[`how`](https://rdrr.io/pkg/permute/man/how.html) for display. The
formatting is more compact than the one used in `print` in the permute
package, and shows only non-default choices. This output is normally
used when printing the results of vegan permutations.

`pasteCall` prints the function call so that it is nicely wrapped in
[`Sweave`](https://rdrr.io/r/utils/Sweave.html) output.

`veganCovEllipse` finds the coordinates for drawing a covariance
ellipse.

`veganMahatrans` transforms data matrix so that its Euclidean distances
are Mahalanobis distances. The input data `x` must be a matrix centred
by columns, and `s2` its covariance matrix. If `s2` is not given,
covariance matrix is found from `x` within the function. If
`na.rm = TRUE`, [`cov`](https://rdrr.io/r/stats/cor.html) is called with
`use = "pairwise.complete.obs"`.

`hierParseFormula` returns a list of one matrix (left hand side) and a
model frame with factors representing hierarchy levels (right hand side)
to be used in
[`adipart`](https://vegandevs.github.io/vegan/reference/adipart.md),
[`multipart`](https://vegandevs.github.io/vegan/reference/multipart.md)
and
[`hiersimu`](https://vegandevs.github.io/vegan/reference/adipart.md).

`GowerDblcen` performs the Gower double centring of a matrix of
dissimilarities. Similar function was earlier available as a compiled
code in stats, but it is not a part of official API, and therefore we
have this poorer replacement.

`addLingoes` and `addCailliez` find the constant added to non-diagonal
(squared) dissimilarities to make all eigenvalues non-negative in
Principal Co-ordinates Analysis
([`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md),
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)).
Function [`cmdscale`](https://rdrr.io/r/stats/cmdscale.html) implements
the Cailliez method. The argument is a matrix of dissimilarities.
