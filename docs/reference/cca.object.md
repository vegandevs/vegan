# Result Object from Constrained Ordination

Ordination methods
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) and
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)
return similar result objects. All these methods use the same internal
function `ordConstrained`. They differ only in (1) initial
transformation of the data and in defining inertia, (2) weighting, and
(3) the use of rectangular rows \\\times\\ columns data or symmetric
rows \\\times\\ rows dissimilarities:
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md) initializes
data to give variance or correlations as inertia,
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) is based on
double-standardized data to give Chi-square inertia and uses row and
column weights,
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md) maps
the real part of dissimilarities to rectangular data and performs RDA,
and [`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md)
performs an RDA-like analysis directly on symmetric dissimilarities.

Function `ordConstrained` returns the same result components for all
these methods, and the calling function may add some more components to
the final result. However, you should not access these result components
directly (using `$`): the internal structure is not regarded as stable
application interface (API), but it can change at any release. If you
access the results components directly, you take a risk of breakage at
any vegan release. The vegan provides a wide set of accessor functions
to those components, and these functions are updated when the result
object changes. This documentation gives an overview of accessor
functions to the `cca` result object.

## Usage

``` r
ordiYbar(x, model = c("CCA", "CA", "pCCA", "partial", "initial"))
# S3 method for class 'cca'
model.frame(formula, ...)
# S3 method for class 'cca'
model.matrix(object, ...)
# S3 method for class 'cca'
weights(object, display = "sites", ...)
```

## Arguments

- object, x, formula:

  A result object from
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
  [`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md), or
  [`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md).

- model:

  Show constrained (`"CCA"`), unconstrained (`"CA"`) or conditioned
  “partial” (`"pCCA"`) results. In `ordiYbar` the value can also be
  `"initial"` for the internal working input data, and `"partial"` for
  the internal working input data after removing the partial effects.

- display:

  Display either `"sites"` or `"species"`.

- ...:

  Other arguments passed to the the function.

## Details

The internal (“working”) form of the dependent (community) data can be
accessed with function `ordiYbar`. The form depends on the ordination
method: for instance, in
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) the data are
weighted and Chi-square transformed, and in
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) they are
Gower-centred dissimilarities. The input data in the original
(“response”) form can be accessed with
[`fitted.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md)
and
[`residuals.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md).
Function
[`predict.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md)
can return either working or response data, and also their lower-rank
approximations.

The model matrix of independent data (“Constraints” and “Conditions”)
can be extracted with `model.matrix`. In partial analysis, the function
returns a list of design matrices called `Conditions` and `Constraints`.
If either component was missing, a single matrix is returned. The
redundant (aliased) terms do not appear in the model matrix. These terms
can be found with
[`alias.cca`](https://vegandevs.github.io/vegan/reference/goodness.cca.md).
Function `model.frame` tries to reconstruct the data frame from which
the model matrices were derived. This is only possible if the original
model was fitted with `formula` and `data` arguments, and still fails if
the `data` are unavailable.

The number of observations can be accessed with
[`nobs.cca`](https://vegandevs.github.io/vegan/reference/nobs.cca.md),
and the residual degrees of freedom with
[`df.residual.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md).
The information on observations with missing values can be accessed with
[`na.action`](https://rdrr.io/r/stats/na.action.html). The terms and
formula of the fitted model can be accessed with
[`formula`](https://rdrr.io/r/stats/formula.html) and
[`terms`](https://rdrr.io/r/stats/terms.html).

The weights used in
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md) can be
accessed with `weights`. In unweighted methods
([`rda`](https://vegandevs.github.io/vegan/reference/cca.md)) all
weights are equal.

The ordination results are saved in separate components for partial
terms, constraints and residual unconstrained ordination. There is no
guarantee that these components will have the same internal names as
currently, and you should be cautious when developing scripts and
functions that directly access these components.

The constrained ordination algorithm is based on QR decomposition of
constraints and conditions (environmental data), and the QR component is
saved separately for partial and constrained components. The QR
decomposition of constraints can be accessed with
[`qr.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md).
This will also include the residual effects of partial terms
(Conditions), and it should be used together with
`ordiYbar(x, "partial")`. The environmental data are first centred in
`rda` or weighted and centred in `cca`. The QR decomposition is used in
many functions that access `cca` results, and it can be used to find
many items that are not directly stored in the object. For examples, see
[`coef.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`coef.rda`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`vif.cca`](https://vegandevs.github.io/vegan/reference/goodness.cca.md),
[`permutest.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`predict.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`predict.rda`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`calibrate.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md).
See [`qr`](https://rdrr.io/r/base/qr.html) for other possible uses of
this component. For instance, the rank of the constraints can be found
from the QR decomposition.

The eigenvalues of the solution can be accessed with
[`eigenvals.cca`](https://vegandevs.github.io/vegan/reference/eigenvals.md).
Eigenvalues are not evaluated for partial component, and they will only
be available for constrained and residual components.

The ordination scores are internally stored as (weighted) orthonormal
scores matrices. These results can be accessed with
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
and
[`scores.rda`](https://vegandevs.github.io/vegan/reference/plot.cca.md)
functions. The ordination scores are scaled when accessed with
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md)
functions, but internal (weighted) orthonormal scores can be accessed by
setting `scaling = FALSE`. Unconstrained residual component has species
and site scores, and constrained component has also fitted site scores
or linear combination scores for sites and biplot scores and centroids
for constraint variables. The biplot scores correspond to the
`model.matrix`, and centroids are calculated for factor variables when
they were used. The scores can be selected by defining the axes, and
there is no direct way of accessing all scores of a certain component.
The number of dimensions can be assessed from
[`eigenvals`](https://vegandevs.github.io/vegan/reference/eigenvals.md).
In addition, some other types can be derived from the results although
not saved in the results. For instance, regression scores and model
coefficients can be accessed with
[`scores`](https://vegandevs.github.io/vegan/reference/scores.md) and
[`coef`](https://rdrr.io/r/stats/coef.html) functions. Partial component
will have no scores.

Distance-based methods
([`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md),
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md)) can
have negative eigenvalues and associated imaginary axis scores. In
addition, species scores are initially missing in
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md) and they
are accessory and found after analysis in
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md) (and
may be misleading). Function
[`sppscores`](https://vegandevs.github.io/vegan/reference/sppscores.md)
can be used to add species scores or replace them with more meaningful
ones.

## See also

The core function is `ordConstrained` which is called by
[`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`rda`](https://vegandevs.github.io/vegan/reference/cca.md),
[`dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.md),
[`capscale`](https://vegandevs.github.io/vegan/reference/dbrda.md) as
well as by unconstrained methods
[`pca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`ca`](https://vegandevs.github.io/vegan/reference/cca.md) and
[`pco`](https://vegandevs.github.io/vegan/reference/dbrda.md). The basic
class is `"cca"` for all methods, and the following functions are
defined for this class:
[`add1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md),
[`alias.cca`](https://vegandevs.github.io/vegan/reference/goodness.cca.md),
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`as.mlm.cca`](https://vegandevs.github.io/vegan/reference/vegan-defunct.md),
[`biplot.cca`](https://vegandevs.github.io/vegan/reference/biplot.rda.md),
[`bstick.cca`](https://vegandevs.github.io/vegan/reference/screeplot.cca.md),
[`calibrate.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`coef.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`cooks.distance.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`deviance.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md),
[`df.residual.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`drop1.cca`](https://vegandevs.github.io/vegan/reference/add1.cca.md),
[`eigenvals.cca`](https://vegandevs.github.io/vegan/reference/eigenvals.md),
[`extractAIC.cca`](https://vegandevs.github.io/vegan/reference/deviance.cca.md),
[`fitted.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`goodness.cca`](https://vegandevs.github.io/vegan/reference/goodness.cca.md),
[`hatvalues.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`labels.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
`model.frame.cca`, `model.matrix.cca`,
[`nobs.cca`](https://vegandevs.github.io/vegan/reference/nobs.cca.md),
[`permutest.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.md),
[`plot.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`points.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`predict.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
`print.cca`,
[`qr.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`residuals.cca`](https://vegandevs.github.io/vegan/reference/predict.cca.md),
[`RsquareAdj.cca`](https://vegandevs.github.io/vegan/reference/RsquareAdj.md),
[`rstandard.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`rstudent.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`screeplot.cca`](https://vegandevs.github.io/vegan/reference/screeplot.cca.md),
[`sigma.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`simulate.cca`](https://vegandevs.github.io/vegan/reference/simulate.rda.md),
[`SSD.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
[`stressplot.cca`](https://vegandevs.github.io/vegan/reference/stressplot.wcmdscale.md),
[`summary.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`text.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.md),
[`tolerance.cca`](https://vegandevs.github.io/vegan/reference/tolerance.md),
[`vcov.cca`](https://vegandevs.github.io/vegan/reference/influence.cca.md),
`weights.cca` . Other functions handling `"cca"` objects include
[`inertcomp`](https://vegandevs.github.io/vegan/reference/goodness.cca.md),
[`intersetcor`](https://vegandevs.github.io/vegan/reference/goodness.cca.md),
[`mso`](https://vegandevs.github.io/vegan/reference/mso.md),
[`ordistep`](https://vegandevs.github.io/vegan/reference/ordistep.md),
[`ordiR2step`](https://vegandevs.github.io/vegan/reference/ordistep.md)
and
[`vif.cca`](https://vegandevs.github.io/vegan/reference/goodness.cca.md).
Functions that can be regarded as special cases of `"cca"` methods
include
[`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.md) and
[`varpart`](https://vegandevs.github.io/vegan/reference/varpart.md).

## Note

The latest large change of result object was made in release 2.5-1 in
2016. You can modernize ancient stray results with
`modernobject <- update(ancientobject)`.

## References

Legendre, P. and Legendre, L. (2012) *Numerical Ecology*. 3rd English
ed. Elsevier.

## Author

Jari Oksanen
