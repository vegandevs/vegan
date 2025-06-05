# vegan 2.7-1

## Installation

* **vegan** no longer depends on **lattice**, but only imports
  **lattice** functions. The **lattice** package is no longer
  automatically loaded. To use **lattice** functions directly, you
  must first attach the package with `library(lattice)`. Longer-term
  plan is to remove **lattice** functions as soon as more modern
  alternatives in **ggplot2** are made available. See Discussion
  [#727](https://github.com/vegandevs/vegan/discussions/727) and
  section Deprecated and Defunct for the changes in this release.

* **vegan** no longer suggests **tcltk**. See `orditkplot` in section
  Deprecated and Defunct.

## New Functions

* Added a set of functions to add new points to an existing NMDS
  ordination from `metaMDS` or `monoMDS`. This serves the same purpose
  as adding new points to an existing eigenvector ordination (for
  instance, `predict.rda`). The main function is `MDSaddpoints`. This
  needs an input of rectangular matrix of dissimilarities of all new
  points (rows) to all old points (columns). Support function
  `dist2xy` can extract needed matrices from dissimilarities of all
  (old and new) points, and function `designdist2` can directly find
  the needed dissimilarities between two data matrices. In addition,
  **analogue** package can calculate such rectangular dissimilarities,
  including many indices that cannot be defined with `designdist2`.

* `betadistances`: new functions to find distances of points to
  group centroids in `betadisper`.  See
  https://stackoverflow.com/questions/77391007/ and issue
  [#606](https://github.com/vegandevs/vegan/issues/606).

* `permulattice`: new function to use lattice graphics for
  `permustats` results without need to first issue `library(lattice)`.

* `optspace`: a new function for matrix completion or filling missing
  elements in a matrix.  The function is used in robust Aitchison
  distance (see below).

## New Features in Ordination Graphics

* `plot.cca` graphics can be configured. `plot.cca` had hard-coded
  graphical parameters and user arguments were ignored (see issue
  [#656](https://github.com/vegandevs/vegan/issues/656)). Now
  graphical parameters can be given either for all score types, or
  with a list of graphical parameters for a specific score.

  The new features are more extensively described in help pages of
  `plot.cca`, `ordiplot` and `biplot.rda`.

* `text.ordiplot` and hence `plot.cca` gained argument `optimize` that
  will call `ordipointlabel` to optimize the location of the text to
  minimize over-writing, but mark the exact score with a point.

* `text.ordiplot` and hence `plot.cca` gained argument `bg=<colour>`
  that will plot text over non-transparent label using `ordilabel`.

* Alternatively ordination plots can be built up adding each score
  type in piped commands. Pipes were available since **vegan** 2.5-1,
  but their use is now improved: `ordilabel` can be used in a pipe,
  `text` can use opaque background label, and `text` and `points`
  function (for `ordiplot`) gained argument for adjusting arrow
  lengths similarly as in `cca`.

* `text.cca` and `points.cca` were completely redesigned because of
  the concerns raised in PR
  [#729](https://github.com/vegandevs/vegan/pull/729). Support
  function `labels.cca` now accepts abbreviated names of score types.

* `text` functions for ordination graphics have arguments `labels` to
  rename `text`and `select` to show only some items. Now these
  functions are consistent and use first `select` and then `labels`
  for the selected subset. Concerns functions `ordilabel`,
  `ordipointlabel`, `orditorp` as well as `text` functions for `cca`
  and friends, `decorana`, `monoMDS`, `metaMDS` and `ordiplot`. See
  issue [#730](https://github.com/vegandevs/vegan/pull/730).

* `orditorp`, `points.cca` and `text.cca` did not accept row names or
  `labels` in `select`. PR
  [#729](https://github.com/vegandevs/vegan/pull/729).

  Species scores can be added to `monoMDS` with `sppscores` function,
  and now these can be accessed in `points` and `text` functions.

* `ordipointlabel` can be used in pipe. Function gained argument
  `label` that allows changing plotted text, and a function `labels`
  that return the current labels. The optimization rules were changed
  to give a slight preference for putting labels outwards from origin
  but avoiding corner positions.
  
* `orditorp` can be used in pipe.

## New Features in Other Functions

* Constrained ordination models (`cca`, `rda`, `dbrda`) inform users
  on completely aliased conditions or constraints, and behave more
  robustly with these degenerate cases. If a model component
  (condition, constrained, residual unconstrained) is completely
  aliased, it still appears in summary table with rank and
  inertia 0. See https://stackoverflow.com/questions/79613784/ and
  issue [#682](https://github.com/vegandevs/vegan/issues/682).

* Robust Aitchison distance uses matrix completion to estimate the
  missing values (`-Inf`) that result from log transformation of zeros
  of the original input data. Earlier we used only above zero values,
  or in simple Aitchison distance replaced zeros with an arbitrary
  pseudocount. For matrix completion **vegan** adds new function
  `optspace` which also can be used independently. The Robust
  Aitchison distance is directly evaluated in `vegdist`, and the
  needed transformation can be performed in `decostand`. PR
  [#667](https://github.com/vegandevs/vegan/pull/667).

* `ordiR2step` calls current model `<model>` instead of `<none>`.

* `vegemite` and `tabasco` can now `use` a factor to show a
  classification. The factor levels and sites within levels can be
  reordered to give a diagonal pattern, as default in `tabasco`
  and in `vegemite` with new argument `diagonalize = TRUE` (defaults
  `FALSE`). With the same argument, `vegemite` can also reorder
  dendrogram (or tree) to give a diagonal pattern. If `coverscale` is
  used, all internal calculations for ordering rows and columns will
  be based on scaled data.

* `make.cepnames` was completely re-designed and is much more flexible
  with enhanced user-control.

* `wascores` can now calculate (unbiased) weighted standard deviation
  of weighted averages with argument `stdev = TRUE`.

* `mantel`, `mrpp`, `oecosimu` and `ordiareatest` gained `summary`
  methods based on `summary.permustats`.

## Bug Fixes

* `ordistep` never dropped aliased terms.

* `permatswap` failed to set some null models (`swsh_samp`,
  `swsh_both`, `backtrack`).

* `df.residual.cca` ignored conditions (partial component). This
  influences many diagnostic statistics documented together with
  `cooks.distance.cca` and `rstandard.cca`.

* `densityplot.permustats` did not know argument `observed` to display
  the observed statistic.

## Deprecated and Defunct

* Disabled use of `summary` to get ordination scores: use `scores` to
  get scores!  For `summary.cca` see
  [#644](https://github.com/vegandevs/vegan/discussions/644).

* lattice function `ordicloud` is deprecated. It is still available in
  CRAN package **vegan3d** (version 1.4-0) as function `ordilattice3d`.

* lattice function `ordisplom` is deprecated: it had bad design and
  was not very useful.

* lattice function `ordiresids` is deprecated: you can access the same
  items using `fitted`, `residuals`, `rstandard`, `rstudent` _etc_ and
  design your own plots.

* `summary.decorana` is defunct. It did nothing useful, but you can
  extract the same information with `scores` and `weights`.

* `orditkplot` was moved to **vegan3d** package and is defunct in
  **vegan**.

* relic function `vegandocs` is officially defunct. Better tools to
  read **vegan** documentation are `browseVignettes("vegan")` and
  `news(package="vegan")`. The function was deprecated in **vegan**
  2\.3-4.

# vegan 2\.6-10

## Startup

* Prints startup message ("This is vegan 2.6-10") only in
  interactive sessions. Version number is no longer shown in
  package checks and other scripts.

## Bug Fixes

* `anova.cca(..., by="margin")` failed when a constraint was
  completely aliased by conditions. See
  [#701](https://github.com/vegandevs/vegan/pull/701).

* `envfit` failed when ordination scores were given in a plain matrix
  instead of a complex ordination result object. Issue
  [#713](https://github.com/vegandevs/vegan/issues/713).

  `envfit` could fail when it was called with only one environmental
  variable *without* formula interface. Formula interface worked
  correctly. Issue [#720](https://github.com/vegandevs/vegan/issues/720).

* `vegemite` dropped dimensions when only one site or species was
  requested.

  `vegemite` could fail with variable lengths of row names (SU
  names).

# vegan 2.6-8

## New Features

* Wrappers for the unconstrained ordination methods principal components
  analysis (PCA), correspondence anslysis (CA), and principal coordinates
  analysis (PCO) are now available via `pca()`, `ca()`, and `pco()`
  respectively. The underlying methods used are `rda()`, `cca()` and `dbrda()`
  respectively. See
  [#655](https://github.com/vegandevs/vegan/issues/655).

* The output from the ordination methods `pca()`, `pco()`, `ca()`,
  `rda()`, `cca()`, `capscale`, and `dbrda()` has changed slightly to
  better separate the results from notifications to the user about
  issues encountered with the data or the model. Related to changes in
  [#682](https://github.com/vegandevs/vegan/issues/682).

* The constrained ordination functions are now louder at informing users when
  one or more terms in a model are aliased (linearly dependent) and their 
  effects cannot be estimated. See
  [#682](https://github.com/vegandevs/vegan/issues/682).

* `cca` and `rda` return centroids for factor levels even when they
  are called without formula, for instance, as `cca(dune, dune.env)`.

* `plot.cca` retains default graphical settings also when only one set
  of scores was displayed.

* `ordiplot` did not pass character size (`cex`) to `plot.cca`. Version
  2.7-0 has more extensive changes, but this fixes the immediate issue
  [#656](https://github.com/vegandevs/vegan/issues/656).

* `adonis2()` now defaults to running an omnibus test of the model
  (`by = NULL`) instead of a sequential test of model terms (`by =
  "terms"`). This makes `adonis2()` more consistent with the default
  for related ordination methods.  See
  [#677](https://github.com/vegandevs/vegan/issues/677).

* `decorana` checks now that input data are numeric instead of
  confusing error message (see
  https://stackoverflow.com/questions/78666646/).

* `make.cepnames` no longer splits names by hyphen: _Capsella
  bursa-pastoris_ used to be `Capspast` but now is `Capsburs`.

## Bug Fixes

* `dbrda` failed in rare cases when an ordination component had only
  negative eigenvalues. Issue
  [#670](https://github.com/vegandevs/vegan/issues/670).

* `plot.cca`: biplot or regression arrows were not nicely scaled and
  drew no arrows when displayed as the only item in graph.

* `ordipointlabel` failed with `decorana` result. Bounding box for
  text could be wrongly estimated with varying values of `cex`.

* `vegdist` with argument `na.rm = TRUE` still failed with missing
  values. Dissimilarity methods `"chisq"` (Chi-square distance) and
  `"mahalanobis"` did not implement `na.rm = TRUE`. Even when missing
  values are removed in calculation, dissimilarities may contain `NA`
  depending on the number and pattern of missing values and
  dissimilarity method.

* `decostand` standardization method `"clr"` did not implement
  `na.rm = TRUE`
  (issue [#661](https://github.com/vegandevs/vegan/issues/661)).
  Standardization methods `"rank"` and `"rrank"` did not retain `NA`
  values but changed them to 0. Original `NA` values are kept in
  `decostand`, but with `na.rm = TRUE` they are ignored when
  transforming other data values.

* `metaMDS`: half-change scaling failed when `maxdist` was fixed, but
  was not 1.

* `summary.ordihull` (and hence `ordiareatest` for convex hulls)
  failed if input had more than two dimensions.

* `simulate.rda` failed with univariate response.

* `vegemite` returned only the last page of multi-page table in its
  (invisible) return object.

# vegan 2.6-6.1

## Bug Fixes

* C function `do_wcentre` (weighted centring) can segfault due to a
  protection error. The problem was found in automatic CRAN
  checks. `do_wcentre` is an internal function that is called from
  `envfit` (`vectorfit`), `wcmdscale` and `varpart` (`simpleCCA`)
  Fixes bug [#653](https://github.com/vegandevs/vegan/issues/653).

# vegan 2.6-6

## Installation

* **vegan** depends on **R** version 4.1.0.

* It is possible to build **vegan** with webR/wasm Fortran
  compiler. Issue [#623](https://github.com/vegandevs/vegan/issues/623).

## New Features

* Permutation tests for CCA were completely redesigned to follow C.J.F
  ter Braak & D.E. te Beest: Environ Ecol Stat 29, 849–868 (2022)
  (https://doi.org/10.1007/s10651-022-00545-4). The constraints are
  now re-weighted for the permuted response data, and in partial model
  they are also residualized by conditions (partial terms). In
  **vegan** (after release 2.4-6) the tests were identical to Canoco,
  but ter Braak & te Beest demonstrated that the results are biased.
  In old **vegan** (release 2.4-2 and earlier) the predictors were
  re-weighted but not residualized. Re-weighting was sufficient to
  remove bias with moderate variation of weights, but residualizing of
  predictors is necessary with strongly varying weights. See
  discussion in issue
  [#542](https://github.com/vegandevs/vegan/issues/542).
  The new scheme only concerns CCA which is a weighted method, and RDA
  and dbRDA permutation is unchanged.

* `summary` of ordination results no longer prints ordination scores
  that often are so voluminous that they hide the real summary; see issue
  [#203](https://github.com/vegandevs/vegan/issues/203). Ordination
  scores should be extracted with `scores` function. This breaks some
  CRAN packages that use `summary.cca` to extract scores. These should
  switch to use `scores`. The maintainers have been contacted and
  patch files are suggested to adapt to this change. See
  [instructions](https://github.com/vegandevs/vegan/discussions/644)
  to fix the packages.

* `scores` function for constrained ordination (CCA, RDA,dbRDA)
  default to return all types of scores (`display = "all"`). Function
  can optionally return a single type of scores as a list of one matrix
  instead of returning a matrix (new argument `droplist`).

* Constrained ordination objects (`cca`, `rda`, `dbrda`) fitted
  without formula interface can have permutation tests (`anova`) by
  `"axis"` and by `"onedf"`. Models by `"terms"` and `"margin"` are only
  possible with formula interface.

* Permutation tests for constrained ordination objects (`cca`, `rda`,
  `dbrda`) with `by = "axis"` stop permutations of later axis once the
  `cutoff` limit is reached. Earlier `cutoff` had to be exceeded. The
  default is to stop permutations once _P_-value 1 is reached. The
  analysis takes care that _P_-values of axes are non-decreasing
  similarly as in Canoco.

* Coefficients of effects in `prc` models are scaled similarly as they
  were scaled in **vegan** pre 2\.5-1. The change was suggested by
  Cajo ter Braak.

* Handling of negative eigenvalues was changed in the `summary` of
  `eigenvals`. Negative eigenvalues are given as negative
  "explanation", and the accumulated proportions add up over 1 for the
  last non-negative eigenvalue, and 1 for the last negative
  eigenvalue.

* The printed output of `capscale` shows proportions for real
  components only and ignores imaginary dimensions. This is consistent
  to `summary` and other support methods. Issue
  [#636](https://github.com/vegandevs/vegan/issues/636).

* `RsquareAdj` of `capscale` is based only on positive eigenvalues,
  and imaginary components are ignored.

* `stressplot.dbrda` refuses to handle partial models. Only the first
  component of variation can be displayed because `dbrda` internal
  ("working") data structures are not additive. For unconstrained
  model `"CA"`, for constrained `"CCA"` and for partial none.

* `predict` for `dbrda` will return the actual
  `type = "working"`. Earlier it returned `"lc"` scores weighted by
  eigenvalues. Both generated same distances and eigenvalues, though.

## Bug Fixes

* Parallel processing was inefficiently implemented and could be
  slower than non-parallel in permutation tests for constrained
  ordination and `adonis2`.

* `plot` and `scores` for `cca` and `rda` family of methods gave an
  error when non-existing axes were requested. Now ignores requests to
  axes numbers that are higher than in the result object.

* `summary` of `prc` ignored extra parameters (such as `const`).

* Over-fitted models with high number of aliased variables caused a
  rare failure in `adonis2` and permutation tests of constrained
  ordination methods (`cca`, `rda`, `dbrda`, `capscale`) with
  arguments `by = "margin"` or `by = "axis"`. This also concerned
  `vif.cca` and `intersetcor`. Typically this occurred with high-order
  interactions of factor variables. See issues
  [#452](https://github.com/vegandevs/vegan/issues/452) and
  [#622](https://github.com/vegandevs/vegan/issues/622)
  
* Some methods accept rectangular raw data input as alternative to
  distances, but did not pass all arguments to distance
  functions. These arguments in `vegdist` could be `binary = TRUE` or
  `pseudocount` with Aitchison distance. This concerns `dbrda`,
  `capscale` and `bioenv`. See issue
  [#631](https://github.com/vegandevs/vegan/issues/631)
  
* `simper` gave arbitrary *p*-values for species that did not occur in
  a subset. Now these are given as `NA`. See
  https://stackoverflow.com/questions/77881877/
  
* `Rsquare.adj` gave arbitrary *p*-values for over-fitted models with
  no residual variation. Now returns `NA` when _R_<sup>2</sup> cannot
  be adjusted. Automatic model building could proceed to such cases,
  and this was fixed in `ordiR2step` which returns _R_<sup>2</sup> = 0
  for overfitted cases. The constrained ordination methods issue a
  warning if the model has no residual component. See issue
  [#610](https://github.com/vegandevs/vegan/issues/610)
  
* `inertcomp(..., display = "sites", proportional = TRUE)` gave wrong
  values.

## Data Sets

* Extended the description of the BCI data sets to avoid
  confusion. The complete BCI survey includes all stems of down to
  1&nbsp;cm DBH, but the BCI data set in **vegan** is a subset of stems of
  DBH 10&nbsp;cm that was published in
  [Science 295, 666&mdash;669, 2002](https://www.science.org/doi/10.1126/science.1066854).
  The data set is intended only to demonstrate methods in **vegan** and for
  ecological research we suggest contacting the BCI team and using the
  complete surveys made available in
  [Dryad](https://doi.org/10.15146/5xcp-0d46).

## Deprecated and Defunct

* `adonis` is deprecated: use `adonis2`. There are several CRAN
  packages that still use `adonis` although we have contacted all
  their authors in June 2022 and again in April 2024, and printed a
  message of forthcoming deprecation since **vegan** 2.6-2. See issue
  [#523](https://github.com/vegandevs/vegan/issues/523). See
  [instructions](https://github.com/vegandevs/vegan/discussions/641)
  to adapt your packages and functions to use `adonis2`.
  
* `orditkplot` was moved to CRAN package **vegan3d** and is deprecated
  in **vegan**. See issue
  [#585](https://github.com/vegandevs/vegan/issues/585) and
  announcement
  [#632](https://github.com/vegandevs/vegan/discussions/632)

* The use of `summary` to extract ordination scores is deprecated: you
  should use `scores` to extract scores. This version still allows
  extracting scores with `summary`, but this will fail in next
  versions. For `summary.cca` see
  [instructions](https://github.com/vegandevs/vegan/discussions/644)
  to change your package.

* Support was removed from ancient `cca` objects (results of `cca`,
  `rda`, `dbrda` or `capscale`) generated before CRAN release 2.5
  (2016). If you still have such stray relics, use
  `newobject <- update(ancientobject)` to modernize the result.

* `as.mcmc.oecosimu` and `as.mcmc.permat` are defunct: use `toCoda`.

* Code of defunct functions was completely removed.

# vegan 2.6-4

## New Features

* Support of `scores` for
  [**ggplot2**](https://CRAN.R-project.org/package=ggplot2) graphics is
  improved and extended for ordination functions. Suitable scores can be
  requested with argument `tidy = TRUE`, and in general all available
  types of scores are returned in a data frame with variable `score`
  labelling the type. The option was implemented in default method of
  `scores` and for structured `wcmdscale` objects, and glitches were
  fixed for `rda` family and `decorana`. Previously `tidy` scores were
  implemented for `cca`, `rda`, `dbrda` family of methods, `metaMDS`,
  `envfit` and `rarecurve`.

* `adonis2` and `anova` for constrained ordination results can perform a
  sequential test of one-degree-of-freedom effects where multi-level
  factors are split to their contrasts. Previously the test was
  available only in `permutest`.

* New `summary` function for `varpart` for a brief overview. The summary
  shows unique and overall contributed variation for each set of
  variables. The fractions shared by several sets of variables are
  divided equally with all contributing sets following Lai J, Zou Y,
  Zhang J, Peres-Neto P (2022\) *Methods in Ecology and Evolution*, 13:
  782–788\.

* `decorana` estimates orthogonalized eigenvalues and the total inertia
  (scaled Chi-square). Orthogonalized eigenvalues can add up to the
  total inertia. Together these enabled implementing `eigenvals`,
  `bstick` and `screeplot` methods for `decorana`.

* Axis lengths are reported for all `decorana` methods.

* Implemented `tolerance` method for `decorana`. This returns the
  criterion that was used in rescaling DCA, and can be used to inspect
  the success of rescaling: it should be constant 1 over the whole axis.

* New `toCoda` function to transform sequential null model results
  from `oecosimu` to an object that can be analysed with
  [**coda**](https://CRAN.R-project.org/package=coda) for convergence
  and independence as an MCMC model. Function replaces
  `as.mcmc.oecosimu` and `as.mcmc.permat`.

* `metaMDS` is more informative about finding similar repeated results
  with random starts and uses less confusing language when reporting the
  results.

* Hellinger distance is directly available in `vegdist`.

* `vegdist`, `betadiver` and `raupcrick` set attribute `maxdist` giving
  the numeric value of theoretical maximum of the dissimilarity index.
  For many dissimilarities this is 1, but &radic;2 for Chord and
  Hellinger distances, for instance. The attribute is `NA` for open
  indices that do not have such a ceiling. `betadiver` has three
  similarity indices and these set `maxdist` 0\.

* `metaMDS` defaults to halfchange scaling when the dissimilarities have
  a numeric `maxdist` attribute, and adapt the threshold to the ceiling
  value. For open indices without ceiling, the threshold will be in the
  scale of dissimilarities. `metaMDS` used a simple test to detect index
  ceiling 1, but the test is now more robust and can also find other
  maximum values. If such inference is made, the function will broadcast
  a message of assumed value of the ceiling.

* Mountford index in `vegdist` is now scaled to maximum value log(2).
  Earlier Mountford distances were scaled to maximum 1\.

* `hatvalues` of constrained ordination objects can sometimes be
  practically 1 or above 1, but now these cases will be exactly 1\. In
  those cases `rstandard`, `rstudent` and `cooks.distance` will be
  `NaN`. The behaviour is similar as in `stats::lm.influence` functions.

* `as.rad` can handle multi-row data frames or matrices and return a
  list of Rank-Abundance data for each row. Earlier only one site was
  handled.

* `decostand` returns attribute `parameters` of settings and variables
  used in standardization. New function `decobackstand` can use
  `parameters` to reconstruct original non-standardized data.
  Back-transformation is not exact but has round-off errors, although
  there is an attempt to keep original zeros exact. Back-transformation
  is not possible for methods `pa`, `rank` and `rrank` and it is not
  implemented for `alr`. Back-transformation queried in
  https://stackoverflow.com/questions/73263526/

* Rarefaction and rarefaction-based methods make sense only with
  original observed counts and give misleading results if data are
  multiplied or rare species are removed. Observed counts usually have
  singletons (species with count one), and these method issue a
  warning if minimum count is higher than one (which may be a false
  positive, but inspect your data). Concerns functions `rarefy`,
  `drarefy`, `rrarefy`, `rarecurve`,
  `specaccum(..., method="rarefy")`, `rareslope` and `avgdist`.
  See github
  [discussion #537](https://github.com/vegandevs/vegan/discussions/537).

* `avgdist` exposes `as.dist` arguments and can return `"dist"`ance
  objects that appear as lower triangles instead of appearing as
  symmetric matrices.

* `betadisper` plots accept `col` argument
  ([PR \#300](https://github.com/vegandevs/vegan/pull/300)).

## Bug Fixes

* `decorana` returned wrong results when Hill's piecewise transformation
  (arguments `before`/`after`) were used, unless downweighting was also
  used.

* `scores` failed when `metaMDS` result had no species scores. Bug was
  introduced in release 2\.6-2\. Issue raised in
  https://stackoverflow.com/questions/72483924/

* `tolerance.cca` failed when only one axis (`choice`) was requested.

* `decostand(..., method="alr")` did not accept name as a `reference`,
  and could fail in some cases.

* CRAN package **proxy** interfered with `simper` and caused an
  obscure error (github issue
  [\#528](https://github.com/vegandevs/vegan/issues/528)).

## Deprecated and Defunct

* `adonis` is on way to deprecation. Use `adonis2` instead.

* `as.mcmc.oecosimu` and `as.mcmc.permat` were deprecated: these could
  not be used as S3 methods without depending on **coda** package. Use
  `toCoda` instead.

# vegan 2.6-2

## Installation

* Compiled code is adapted to the changes in **R** 4\.2\.0\. See issues
  [\#447](https://github.com/vegandevs/vegan/issues/447),
  [\#507](https://github.com/vegandevs/vegan/issues/507).

* Cross-references to function in other packages were adapted to more
  stringent tests in CRAN

## New Features

* Aitchison and robust Aitchison distances were added to `vegdist`.
  Similar data transformations were added to `decostand`.

* Several functions can return “tidy” data structures that can be used
  in ggplot2 graphics: `rarecurve`, `scores` functions for constrained
  ordination (`cca` etc.), `decorana`, `envfit`, `metaMDS`.

* `scores.envfit` gained argument `arrow.mul`. vegan `plot` functions
  used this automatically, but now it is easier to use `envfit` in
  non-vegan plotting.

* Added function `simpson.unb` for unbiased Simpson diversity that is
  more robust to the variation in sample sizes.

* `diversity` gained argument `group` to calculate indices for pooled
  data. Discussed in issue
  [\#393](https://github.com/vegandevs/vegan/issues/393).

* `simper` is much faster even though parallel processing is not
  implemented in the new code.

* `pairs` function was added to plot `permustats` variables against each
  other.

* `varpart` accepts dissimilarities given as a symmetric square matrix
  instead of `"dist"` object per wish of issue
  [\#497](https://github.com/vegandevs/vegan/issues/497).

* `metaMDS` adopted a more user-friendly policy, and `trymax` will
  always be the maximum number of tries. See dicussion in
  https://stackoverflow.com/questions/66748605/.

* `adonis2` accepts `strata`. `adonis2` is the new main function that
  replaces old `adonis`. See issue
  [\#427](https://github.com/vegandevs/vegan/issues/427).

* Fisher alpha (`fisherfit`) is badly suited for extreme communities
  that do not follow Fisher's model. Now `fisherfit` returns `NA` to
  communities that have 0 or 1 species, and issues a warning with
  communities consisting of singletons and having extreme Fisher alpha.

* `adipart` and `multipart` formulae will automatically add unique id
  and and constant. This will always sandwich the requested grouping
  between alpha and gamma diversities, but not change the results for
  requested groupings.

## Bug Fixes

* `anova` function failed in marginal tests when constrained partial
  ordination model (`cca`, `rda` etc.) had interaction terms. Issue
  [\#463](https://github.com/vegandevs/vegan/issues/463).

* Constrained ordination (`cca` etc.) gave misleading results when all
  external variables (constraints, condition) were constant and
  explained nothing.

* `decorana` could fail when some axes had zero eigenvalues. Issue
  [\#401](https://github.com/vegandevs/vegan/issues/401).

* Species accumulation (`specaccum`) failed when there was only one
  species, but several “communities”. Issue
  [\#501](https://github.com/vegandevs/vegan/issues/501).

* Parallel processing failed in Windows or with socket clusters in
  `permutest` of `betadisper`. Issue
  [\#369](https://github.com/vegandevs/vegan/issues/369).

* `orditorp` failed if numeric labels were supplied. Reported in
  https://stackoverflow.com/questions/69272366/.

* Argument `summarize` was accidentally dropped from `goodness.cca` in
  2017\.

* `taxa2dist` failed if there was only one usable taxonomic level. See
  https://stackoverflow.com/questions/67231431/.

## Deprecated and Defunct

* Function `adonis2` will replace `adonis`.

* `humpfit` functions are defunct and removed. They are available in
  non-CRAN package natto at <https://github.com/jarioksa/natto>.

* `commisimulator` is defunct. Use `simulate` for `nullmodel` objects.

* `permuted.index` is finally defunct (it was deprecated in vegan
  2\.2-0\).

* `as.mlm` is defunct. Use functions documented with `influence.cca`,
  such as `hatvalues.cca`, `rstandard.cca`, `rstudent.cca`,
  `cooks.distance.cca` and others.

# vegan 2.5-7

## Bug Fixes

* Several distance-based functions failed if all distances were zero
  (`betadisper`, `capscale`, `isomap`, `monoMDS`, `pcnm`, `wcmdscale`).
  Reported in github issue
  [\#372](https://github.com/vegandevs/vegan/issues/372).

* Non-linear self-starting regression models `SSarrhenius`, `SSgitay`,
  `SSgleason` and `SSlomolino` failed in future **R**. The failure was
  caused by internal changes in **R**-devel. Github issue
  [\#382](https://github.com/vegandevs/vegan/issues/382).

* Arrow labels were in wrong position in `plot.envfit(..., add =
  FALSE)`.

* `rarecurve` added unnecessary names to the results. Github issue
  [\#352](https://github.com/vegandevs/vegan/issues/352).

* `permutest` for `betadisper` failed in parallel processing in Windows
  and in other systems when socket clusters were used. Github issue
  [\#369](https://github.com/vegandevs/vegan/issues/369).

## New Features

* Chi-square and Chord distances were added to `vegdist`. Both of these
  distances can be calculated as Euclidean distances of transformed
  data, and actually were available earlier, but many users did not
  notice this.

* `monoMDS` (and hence `metaMDS`) uses stricter convergence criteria.
  This improves possibilities to find stable solutions. However, users
  may still need to tweak convergence criteria with their data. See
  discussion in Github issue
  [\#354](https://github.com/vegandevs/vegan/issues/354).

* `text` functions for constrained ordination plots (`cca`, `rda`,
  `dbrda`, `capscale`) accept now expression labels. This allows using
  subscripts, superscripts and mathematical expressions. New support
  function `labels.cca` returns the current text labels so that authors
  can change the desired ones. See github issue
  [\#374](https://github.com/vegandevs/vegan/issues/374).

* `vegemite` returns invisibly the final formatted table allowing
  further processing.

* `ordiplot` passes `cex` argument to `linestack` and `decorana` plots.

# vegan 2.5-6

## Bug Fixes

* `vegdist` silently accepted missing values (`NA`) and removed them
  from the analysis also with option `na.rm = FALSE`. The behaviour was
  introduced in vegan version 2\.5-1\. See GitHub issue
  [\#319](https://github.com/vegandevs/vegan/issues/319).

* The labels were displaced when the bunch of arrows was not drawn at
  the origin of the ordination graph in `envfit`. See GitHub issue
  [\#315](https://github.com/vegandevs/vegan/issues/315).

* Hill scale in `coverscale` is open-ended and is not limited to percent
  data, unlike most traditional cover class scales which are undefined
  above 100% cover.

## New Features

* The results of `as.rad` no longer print the index attribute: the
  attribute is still in the object, but printing made the output messy.

# vegan 2.5-5

## Installation

* vegan depends on **R** 3\.4\.0 or higher. The next vegan release may
  increase the dependence to **R** 3\.6\.0\.

* **R** 3\.6\.0 improved the method to find random indices for permuting
  and sampling data. Vegan relies now on the **R** functions in its
  ecological null models (functions `nullmodel`, `oecosimu`, `commsim`,
  `permatfull`, `permatswap` and others). Technically this change is
  compatible with **R** 3\.4\.0 and later, but you can only gain the
  benefits of improved code with a current release of **R**. The null
  models may change due to this change, and most certainly they change
  in **R** 3\.6\.0\. See NEWS for the **R** 3\.6\.0 release and
  discussion in github issue
  [\#312](https://github.com/vegandevs/vegan/issues/312).
  
  Most vegan permutation routines rely on
  [permute](https://CRAN.R-project.org/package=permute), and there you
  gain similar benefits of improved randomness when you upgrade **R**.

* Thanks to the new **R** dependence, `sigma` for constrained ordination
  results works without workarounds of vegan 2\.5-2\. This fixes
  completely the issue discussed in
  [\#274](https://github.com/vegandevs/vegan/issues/274).

* Vegan test results cannot be reproduced in older versions than **R**
  3\.6\.0\. If you are worried about this, you should upgrade **R**.

## Bug Fixes

* `metaMDS` failed in scaling results when other `engine` than `monoMDS`
  was used. However, we recommend you use `monoMDS`. See github issue
  [\#310](https://github.com/vegandevs/vegan/issues/310).

## New Features

* `betadisper` changed interpretation of negative squared distances
  which give complex-valued distances. Now they are regarded as
  zero-distances whereas earlier we used their modulus. This will change
  the results in cases where you had negative squared distances. For
  further discussion, see github issue
  [\#306](https://github.com/vegandevs/vegan/issues/306).

# vegan 2.5-4

## Installation and Testing

* The code for interpreting formula will change in **R** 3\.6\.0, and
  this makes constrained ordination methods (`cca`, `rda`, `dbrda`,
  `capscale`) to fail. See github issue
  [\#299](https://github.com/vegandevs/vegan/issues/299).

* **R** 3\.6\.0 introduces a new environment variable
  `_R_CHECK_LENGTH_1_LOGIC2_`, and several functions fail if this
  variable is set. Changes concern `ordiplot`, `plot` and `summary` for
  constrained ordination objects, and `ordixyplot`. See github issue
  [\#305](https://github.com/vegandevs/vegan/issues/305).

## Bug Fixes

* `decorana` gave incorrect results when downweighting was used
  (argument `iweigh = 1`). The bug was introduced in vegan 2\.5-1 and
  reported as github issue
  [\#303](https://github.com/vegandevs/vegan/issues/303).

* `goodness` for constrained ordination methods failed when the
  constraints had rank = 1 (only one constraining variable). Reported by
  Pierre Legendre.

## New Features

* Adjusted _R_<sup>2</sup> is enabled for partial RDA models (functions `rda` and
  `dbrda`) and partial CCA models (function `cca`) in function
  `RsquareAdj`. The feature was disabled in vegan 2\.5-1 for both. For
  RDA, the calculation is similar as in vegan 2\.4-6 and earlier.
  Partial CCA is now consistent with RDA and differs from the earlier
  implementation. For both methods, the partial models are consistent
  with `varpart`. See github issue
  [\#295](https://github.com/vegandevs/vegan/issues/295).

# vegan 2.5-3

## Installation

* Tests for numerical analysis were written more robustly so that they
  give more similar results with alternative platforms and versions of
  **R** and BLAS/Lapack libraries. See github issue
  [\#282](https://github.com/vegandevs/vegan/issues/282).

## Bug Fixes

* Constrained ordination gave misleading results when some constraints
  or conditions had data with NULL variables. This rarely happens in
  normal usage, but could happen in marginal `anova` as reported in
  github issue [\#291](https://github.com/vegandevs/vegan/issues/291).

* Several functions for numerical analysis wrongly accepted non-numeric
  data (for instance, factors) and gave either meaningless results or
  confusing error messages. Fixed functions include `beals`,
  `designdist`, `diversity`, `gdispweight`, `indpower`, `spantree`,
  `specpool`, `tsallis`, `tsallisaccum` and `vegdist`. See github issue
  [\#292](https://github.com/vegandevs/vegan/issues/292).

* `envfit` with vectors could fail with missing data.

* The original data were not scaled and centred similarly as simulations
  in `simulate.rda` when several simulations were returned as a `simmat`
  object (which is compatible with `nullmodel` simulations and can be
  used in `oecosimu`).

## New Features

* `anosim` checks its input to avoid confusing error messages like that
  reported in https://stackoverflow.com/questions/52082743/.

* Broken-stick distribution (function `bstick`) is no longer calculated
  for distance-based Redundancy Analysis (`dbrda`) with negative
  eigenvalues, because it is not clear how this should be done. Now
  `dbrda` and `capscale` are similar with this respect.

* `print` function for `betadisper` results gained new argument `neigen`
  to select the number of eigenvalues shown. The `print` is more robust
  when the number of eigenvalues is lower than the requested `neigen`.

## Deprecated

* Function `humpfit` was moved to the natto package and is still
  available from <https://github.com/jarioksa/natto>. It is scheduled
  for complete removal in vegan 2\.6-0\.

# vegan 2.5-2

## Installation and Compatibility

* Vegan declares dependence on **R** version 3\.2\.0\. This dependence
  was not yet noticed in the previous vegan release. However, the
  generic `sigma` function was only defined in **R**-3\.3\.0, and
  therefore `sigma.cca` of vegan must be spelt out completely when using
  **R**-3\.2\.x. See discussion in issue
  [\#274](https://github.com/vegandevs/vegan/issues/274).

* CRAN package [klaR](https://CRAN.R-project.org/package=klaR) has
  function `rda`, and when loaded together with vegan this clashes with
  vegan `rda` for Redundancy Analysis. Vegan tries to mitigate the
  problem. In most cases vegan functions will be used if vegan was
  loaded after klaR, and an error message is issued if klaR objects are
  handled with vegan functions. klaR is also tricked to print an
  informative message if it handles vegan objects. However, vegan
  namespace can be attached automatically at the start-up and then klaR
  functions will take precedence. This was reported as issue
  [\#277](https://github.com/vegandevs/vegan/issues/277).

* Bioconductor package phyloseq has a problem with `vegdist` function
  for dissimilarities. The problem can be fixed by re-installing
  phyloseq from its *source package*. If you cannot do this, you must
  either downgrade to vegan version 2\.4-6 or wait till Bioconductor
  binary packages are upgraded. This was reported in
  https://stackoverflow.com/questions/49882886/, and as
  vegan issue [\#272](https://github.com/vegandevs/vegan/issues/272),
  and as phyloseq issues
  [\#918](https://github.com/joey711/phyloseq/issues/918) and
  [\#921](https://github.com/joey711/phyloseq/issues/921).

## Bug Fixes

* Plotting `betadisper` failed if any of the `groups` had only one
  member. Reported in https://stackoverflow.com/questions/50267430/.

* Permutation tests for constrained ordination (`anova.cca`,
  `permutest.cca`) could fail in parallel processing with socket
  clusters. Socket clusters are always used in Windows and they can also
  be used in other operating systems when created with `makeCluster`.
  See issue [\#276](https://github.com/vegandevs/vegan/issues/276).

# vegan 2.5-1

## GENERAL

* This is a major new release with changes all over the package: Nearly
  40% of program files were changed from the previous release. Please
  report regressions and other issues in
  <https://github.com/vegandevs/vegan/issues/>.

* Compiled code is used much more extensively, and most compiled
  functions use `.Call` interface. This gives smaller memory footprint
  and is also faster. In wall clock time, the greatest gains are in
  permutation tests for constrained ordination methods (`anova.cca`) and
  binary null models (`nullmodel`).

* Constrained ordination functions (`cca`, `rda`, `dbrda`, `capscale`)
  are completely rewritten and share most of their code. This makes them
  more consistent with each other and more robust. The internal
  structure changed in constrained ordination objects, and scripts may
  fail if they try to access the result object directly. There never was
  a guarantee for unchanged internal structure, and such scripts should
  be changed and they should use the provided support functions to
  access the result object (see documentation of `cca.object` and github
  issue [\#262](https://github.com/vegandevs/vegan/issues/262)). Some
  support and analysis functions may no longer work with result objects
  created in previous vegan versions. You should use
  `update(old.result.object)` to fix these old result objects. See
  github issues [\#218](https://github.com/vegandevs/vegan/issues/218),
  [\#227](https://github.com/vegandevs/vegan/issues/227).

* vegan includes some tests that are run when checking the package
  installation. See github issues
  [\#181](https://github.com/vegandevs/vegan/issues/181),
  [\#271](https://github.com/vegandevs/vegan/issues/271).

* The informative messages (warnings, notes and error messages) are
  cleaned and unified which also makes possible to provide translations.

## New Functions

* `avgdist`: new function to find averaged dissimilarities from several
  random rarefactions of communities. Code by Geoffrey Hannigan. See
  github issues [\#242](https://github.com/vegandevs/vegan/issues/242),
  [\#243](https://github.com/vegandevs/vegan/issues/243),
  [\#246](https://github.com/vegandevs/vegan/issues/246).

* `chaodist`: new function that is similar to `designdist`, but uses
  Chao terms that are supposed to take into account the effects of
  unseen species (Chao et al., *Ecology Letters* **8,** 148-159; 2005\).
  Earlier we had Jaccard-type Chao dissimilarity in `vegdist`, but the
  new code allows defining any kind of Chao dissimilarity.

* New functions to find influence statistics of constrained ordination
  objects: `hatvalues`, `sigma`, `rstandard`, `rstudent`,
  `cooks.distance`, `SSD`, `vcov`, `df.residual`. Some of these could be
  earlier found via `as.mlm` function which is deprecated. See github
  issue [\#234](https://github.com/vegandevs/vegan/issues/234).

* `boxplot` was added for `permustats` results to display the
  (standardized) effect sizes.

* `sppscores`: new function to add or replace species scores in
  distance-based ordination such as `dbrda`, `capscale` and `metaMDS`.
  Earlier `dbrda` did not have species scores, and species scores in
  `capscale` and `metaMDS` were based on raw input data which may not be
  consistent with the used dissimilarity measure. See github issue
  [\#254](https://github.com/vegandevs/vegan/issues/254).

* `cutreeord`: new function that is similar to `stats::cutree`, but
  numbers the cluster in the order they appear in the dendrogram (left
  to right) instead of labelling them in the order they appeared in the
  data.

* `sipoo.map`: a new data set of locations and sizes of the islands in
  the Sipoo archipelago bird data set `sipoo`.

## New Features in Constrained Ordination

* The inertia of Correspondence Analysis (`cca`) is called “scaled
  Chi-square” instead of using a name of a little known statistic.

* If elements for Constraints and Conditions are data frames in
  non-formula call of `rda` or `cca`, these are automatically expanded
  to model matrices and can contain factor variables. Earlier they had
  to be numerical model matrices and factors could only be used with the
  formula interface.

* Regression scores for constraints can be extracted and plotted for
  constrained ordination methods. See github issue
  [\#226](https://github.com/vegandevs/vegan/issues/226).

* Full model (`model = "full"`) is again enabled in permutations tests
  for constrained ordination results in `anova.cca` and `permutest.cca`.

* `permutest.cca` gained a new option `by = "onedf"` to perform tests by
  sequential one degree-of-freedom contrasts of factors. This option is
  not (yet) enabled in `anova.cca`.

* The permutation tests are more robust, and most scoping issues should
  have been fixed.

* Permutation tests use compiled C code and they are much faster. See
  github issue [\#211](https://github.com/vegandevs/vegan/issues/211).

* `permutest` printed layout is similar to `anova.cca`.

* `eigenvals` gained a new argument (`model`) to select either
  constrained or unconstrained scores. The old argument (`constrained`)
  is deprecated. See github issue
  [\#207](https://github.com/vegandevs/vegan/issues/207).

* `summary.eigenvals` returns a matrix instead of a list containing only
  that matrix.

* Adjusted _R_<sup>2</sup> is not calculated for partial ordination, because it is
  unclear how this should be done (function `RsquareAdj`).

* `ordiresids` can display standardized and studentized residuals.

* Function to construct `model.frame` and `model.matrix` for constrained
  ordination are more robust and fail in fewer cases.

* `goodness` and `inertcomp` for constrained ordination result object no
  longer has an option to find distances: only explained variation is
  available.

* `inertcomp` gained argument `unity`. This will give “local
  contributions to beta-diversity” (LCBD) and “species contribution to
  beta-diversity” (SCBD) of Legendre & De Cáceres (*Ecology Letters*
  **16,** 951-963; 2012\).

* `goodness` is disabled for `capscale`.

* `prc` gained argument `const` for general scaling of results similarly
  as in `rda`.

* `prc` uses regression scores for Canoco-compatibility.

## New Features in Null Model Communities

* The C code for swap-based binary null models was made more efficients,
  and the models are all faster. Many of these models selected a 2
  times 2 submatrix, and for this they generated four random numbers
  (two rows, two columns). Now we skip selecting third or fourth random
  number if it is obvious that the matrix cannot be swapped. Since most
  of time was used in generating random numbers in these functions, and
  most candidates were rejected, this speeds up functions. However, this
  also means that random number sequences change from previous vegan
  versions, and old binary model results cannot be replicated exactly.
  See github issues
  [\#197](https://github.com/vegandevs/vegan/issues/197),
  [\#255](https://github.com/vegandevs/vegan/issues/255) for details and
  timing.

* Ecological null models (`nullmodel`, `simulate`, `make.commsim`,
  `oecosimu`) gained new null model `"greedyqswap"` which can radically
  speed up quasi-swap models with minimal risk of introducing bias.

* Backtracking is written in C and it is much faster. However,
  backtracking models are biased, and they are provided only because
  they are classic legacy models.

## New Features in Other Functions

* `adonis2` gained a column of _R_<sup>2</sup> similarly as old `adonis`.

* Great part of **R** code for `decorana` is written in C which makes it
  faster and reduces the memory footprint.

* `metaMDS` results gained new `points` and `text` methods.

* `ordiplot` and other ordination `plot` functions can be chained with
  their `points` and `text` functions allowing the use of
  [magrittr](https://CRAN.R-project.org/package=magrittr) pipes. The
  `points` and `text` functions gained argument to draw arrows allowing
  their use in drawing biplots or adding vectors of environmental
  variables with `ordiplot`. Since many ordination `plot` methods return
  an invisible `"ordiplot"` object, these `points` and `text` methods
  also work with them. See github issue
  [\#257](https://github.com/vegandevs/vegan/issues/257).

* Lattice graphics (`ordixyplot`) for ordination can add polygons that
  enclose all points in the panel and complete data.

* `ordicluster` gained option to suppress drawing in plots so that it
  can be more easily embedded in other functions for calculations.

* `as.rad` returns the index of included taxa as an attribute.

* Random rarefaction (function `rrarefy`) uses compiled C code and is
  much faster.

* `plot` of `specaccum` can draw short horizontal bars to vertical error
  bars. See https://stackoverflow.com/questions/45378751.

* `decostand` gained new standardization methods `rank` and `rrank`
  which replace abundance values by their ranks or relative ranks. See
  github issue [\#225](https://github.com/vegandevs/vegan/issues/225).

* Clark dissimilarity was added to `vegdist` (this cannot be calculated
  with `designdist`).

* `designdist` evaluates minimum terms in compiled code, and the
  function is faster than `vegdist` also for dissimilarities using
  minimum terms. Although `designdist` is usually faster than `vegdist`,
  it is numerically less stable, in particular with large data sets.

* `swan` passes `type` argument to `beals`.

* `tabasco` can use traditional cover scale values from function
  `coverscale`. Function `coverscale` can return scaled values as
  integers for numerical analysis instead of returning characters for
  printing.

* `varpart` can partition Chi-squared inertia of correspondence
  analysis with new argument `chisquare`. The adjusted _R_<sup>2</sup>
  is based on permutation tests, and the replicate analysis will have
  random variation.

* The explanatory tables can be data frames with factors or single
  factors in `varpart` and these will be automatically expanded to model
  matrices. Earlier factors could only be used with one-sided model
  formulae. Based on the code suggested by Daniel Borcard, Univ.
  Montréal.

## Bug Fixes

* Very long `Condition()` statements (\> 500 characters) failed in
  partial constrained ordination models (`cca`, `rda`, `dbrda`,
  `capscale`). The problem was detected in
  https://stackoverflow.com/questions/49249816.

* Labels were not adjusted when arrows were rescaled in `envfit` plots.
  See https://stackoverflow.com/questions/49259747.

* `ordiArrowMul` failed if there was only one arrow to be plotted in
  `envfit`.

## Deprecated and Defunct

* `as.mlm` function for constrained correspondence analysis is
  deprecated in favour of new functions that directly give the influence
  statistics. See github issue
  [\#234](https://github.com/vegandevs/vegan/issues/234).

* `commsimulator` is now defunct: use `simulate` for `nullmodel`
  objects.

* [ade4](https://CRAN.R-project.org/package=ade4) `cca` objects are no
  longer handled in vegan: ade4 has had no `cca` since version 1\.7-8
  (August 9, 2017\).

# vegan 2.4-6

## Installation and Building

* CRAN packages are no longer allowed to use FORTRAN input, but
  `read.cep` function used FORTRAN format to read legacy CEP and Canoco
  files. To avoid NOTEs and WARNINGs, the function was re-written in
  **R**. The new `read.cep` is less powerful and more fragile, and can
  only read data in “condensed” format, and it can fail in several cases
  that were successful with the old code. The old FORTRAN-based function
  is still available in
  [cepreader](https://CRAN.R-project.org/package=cepreader). See github
  issue [\#263](https://github.com/vegandevs/vegan/issues/263). The
  cepreader package is developed in
  <https://github.com/vegandevs/cepreader>.

## Bug Fixes

* Some functions for rarefaction (`rrarefy`), species abundance
  distribution (`preston`) and species pool (`estimateR`) need exact
  integer data, but the test allowed small fuzz. The functions worked
  correctly with original data, but if data were transformed and then
  back-transformed, they would pass the integer test with fuzz and give
  wrong results. For instance, `sqrt(3)^2` would pass the test as 3, but
  was interpreted strictly as integer 2\. See github issue
  [\#259](https://github.com/vegandevs/vegan/issues/259).

## New Features

* `ordiresids` uses now weighted residuals for `cca` results.

# vegan 2.4-5

## Bug Fixes

* Several “Swap & Shuffle” null models generated wrong number of initial
  matrices. Usually they generated too many, which was not dangerous,
  but it was slow. However, random sequences will change with this fix.

* Lattice graphics for ordination (`ordixyplot` and friends) colour the
  arrows by `groups` instead of randomly mixed colours.

* Information on constant or mirrored permutations was omitted when
  reporting permutation tests (e.g., in `anova` for constrained
  ordination).

## New Features

* `ordistep` has improved interpretation of `scope`: if the lower scope
  is missing, the formula of the starting solution is taken as the lower
  scope instead of using an empty model. See
  https://stackoverflow.com/questions/46985029/.

* `fitspecaccum` gained new support functions `nobs` and `logLik` which
  allow better co-operation with other packages and functions. See
  GitHub issue [\#250](https://github.com/vegandevs/vegan/issues/250).

* The “backtracking” null model for community simulation is faster.
  However, “backtracking” is a biased legacy model that should not be
  used except in comparative studies.

# vegan 2.4-4

## Installation and Building

* `orditkplot` should no longer give warnings in CRAN tests.

## Bug Fixes

* `anova(..., by = "axis")` for constrained ordination (`cca`, `rda`,
  `dbrda`) ignored partial terms in `Condition()`.

* `inertcomp` and `summary.cca` failed if the constrained component was
  defined, but explained nothing and had zero rank. See
  https://stackoverflow.com/questions/43683699/.

* Labels are no longer cropped in the `meandist` plots.

## New Features

* The significance tests for the axes of constrained ordination use now
  forward testing strategy. More extensive analysis indicated that the
  previous marginal tests were biased. This is in conflict with
  Legendre, Oksanen & ter Braak, *Methods Ecol Evol* **2,** 269–277
  (2011\) who regarded marginal tests as unbiased.

* Canberra distance in `vegdist` can now handle negative input entries
  similarly as latest versions of **R**.

# vegan 2.4-3

## Installation and Building

* vegan registers native **C** and **Fortran** routines. This avoids
  warnings in model checking, and may also give a small gain in speed.

* Future versions of vegan will deprecate and remove elements
  `pCCA$Fit`, `CCA$Xbar`, and `CA$Xbar` from `cca` result objects. This
  release provides a new function `ordiYbar` which is able to construct
  these elements both from the current and future releases. Scripts and
  functions directly accessing these elements should switch to
  `ordiYbar` for smooth transition.

## Bug Fixes

* `as.mlm` methods for constrained ordination include zero intercept to
  give the correct residual degrees of freedom for derived statistics.

* `biplot` method for `rda` passes `correlation` argument to the scaling
  algorithm.

* Biplot scores were wrongly centred in `cca` which caused a small error
  in their values.

* Weighting and centring were corrected in `intersetcor` and `spenvcor`.
  The fix can make a small difference when analysing `cca` results.
  
  Partial models were not correctly handled in `intersetcor`.

* `envfit` and `ordisurf` functions failed when applied to species
  scores.

* Non-standard variable names can be used within `Condition()` in
  partial ordination. Partial models are used internally within several
  functions, and a problem was reported by Albin Meyer (Univ Lorraine,
  Metz, France) in `ordiR2step` when using a variable name that
  contained a hyphen (which was wrongly interpreted as a minus sign in
  partial ordination).

* `ordispider` did not pass graphical arguments when used to show the
  difference of LC and WA scores in constrained ordination.

* `ordiR2step` uses only `forward` selection to avoid several problems
  in model evaluation.

* `tolerance` function could return `NaN` in some cases when it should
  have returned `0`. Partial models were not correctly analysed.
  Misleading (non-zero) tolerances were sometimes given for species that
  occurred only once or sampling units that had only one species.

# vegan 2.4-2

## Bug Fixes

* Permutation tests (`permutests`, `anova`) for the first axis failed in
  constrained distance-based ordination (`dbrda`, `capscale`). Now
  `capscale` will also throw away negative eigenvalues when first
  eigenvalues are tested. All permutation tests for the first axis are
  now faster. The problem was reported by Cleo Tebby and the fixes are
  discussed in GitHub issue
  [\#198](https://github.com/vegandevs/vegan/issues/198) and pull
  request [\#199](https://github.com/vegandevs/vegan/pull/199).

* Some support functions for `dbrda` or `capscale` gave results or some
  of their components in wrong scale. Fixes in `stressplot`, `simulate`,
  `predict` and `fitted` functions.

* `intersetcor` did not use correct weighting for `cca` and the results
  were slightly off.

* `anova` and `permutest` failed when `betadisper` was fitted with
  argument `bias.adjust = TRUE`. Fixes Github issue
  [\#219](https://github.com/vegandevs/vegan/issues/219) reported by
  Ross Cunning, O'ahu, Hawaii.

* `ordicluster` should return invisibly only the coordinates of internal
  points (where clusters or points are joined), but last rows contained
  coordinates of external points (ordination scores of points).

* The `cca` method of `tolerance` was returning incorrect values for all
  but the second axis for sample heterogeneities and species tolerances.
  See issue [\#216](https://github.com/vegandevs/vegan/issues/216) for
  details.

## New Features

* Biplot scores are scaled similarly as site scores in constrained
  ordination methods `cca`, `rda`, `capscale` and `dbrda`. Earlier they
  were unscaled (or more technically, had equal scaling on all axes).

* `tabasco` adds argument to `scale` the colours by rows or columns in
  addition to the old equal scale over the whole plot. New arguments
  `labRow` and `labCex` can be used to change the column or row labels.
  Function also takes care that only above-zero observations are
  coloured: earlier tiny observed values were merged to zeros and were
  not distinct in the plots.

* Sequential null models are somewhat faster (up to 10%). Non-sequential
  null models may be marginally faster. These null models are generated
  by function `nullmodel` and also used in `oecosimu`.

* `vegdist` is much faster. It used to be clearly slower than
  `stats::dist`, but now it is nearly equally fast for the same
  dissimilarity measure.

* Handling of `data=` in formula interface is more robust, and messages
  on user errors are improved. This fixes points raised in Github issue
  [\#200](https://github.com/vegandevs/vegan/issues/200).

* The families and orders in `dune.taxon` were updated to APG IV (*Bot J
  Linnean Soc* **181,** 1–20; 2016\) and a corresponding classification
  for higher levels (Chase & Reveal, *Bot J Linnean Soc* **161,**
  122-127; 2009\).

# vegan 2.4-1

## Installation

* Fortran code was modernized to avoid warnings in latest **R**. The
  modernization should have no visible effect in functions. Please
  report all suspect cases as [vegan
  issues](https://github.com/vegandevs/vegan/issues/).

## Bug Fixes

* Several support functions for ordination methods failed if the
  solution had only one ordination axis, for instance, if there was only
  one constraining variable in CCA, RDA and friends. This concerned
  `goodness` for constrained ordination, `inertcomp`, `fitted` for
  `capscale`, `stressplot` for RDA, CCA (GitHub issue
  [\#189](https://github.com/vegandevs/vegan/issues/189)).

* `goodness` for CCA & friends ignored `choices` argument (GitHub issue
  [\#190](https://github.com/vegandevs/vegan/issues/190)).

* `goodness` function did not consider negative eigenvalues of db-RDA
  (function `dbrda`).

* Function `meandist` failed in some cases when one of the groups had
  only one observation.

* `linestack` could not handle expressions in `labels`. This regression
  is discussed in GitHub issue
  [\#195](https://github.com/vegandevs/vegan/issues/195).

* Nestedness measures `nestedbetajac` and `nestedbetasor` expecting
  binary data did not cope with quantitative input in evaluating
  Baselga's matrix-wide Jaccard or Sørensen dissimilarity indices.

* Function `as.mcmc` to cast `oecosimu` result to an MCMC object
  ([coda](https://CRAN.R-project.org/package=coda) package) failed if
  there was only one chain.

## New Features

* `diversity` function returns now `NA` if the observation had `NA`
  values instead of returning `0`. The function also checks the input
  and refuses to handle data with negative values. GitHub issue
  [\#187](https://github.com/vegandevs/vegan/issues/187).

* `rarefy` function will work more robustly in marginal case when the
  user asks for only one individual which can only be one species with
  zero variance.

* Several functions are more robust if their factor arguments contain
  missing values (`NA`): `betadisper`, `adipart`, `multipart`,
  `hiersimu`, `envfit` and constrained ordination methods `cca`, `rda`,
  `capscale` and `dbrda`. GitHub issues
  [\#192](https://github.com/vegandevs/vegan/issues/192) and
  [\#193](https://github.com/vegandevs/vegan/issues/193).

# vegan 2.4-0

## Distance-based Analysis

* Distance-based methods were redesigned and made consistent for
  ordination (`capscale`, new `dbrda`), permutational ANOVA (`adonis`,
  new `adonis2`), multivariate dispersion (`betadisper`) and variation
  partitioning (`varpart`). These methods can produce negative
  eigenvalues with several popular semimetric dissimilarity indices, and
  they were not handled similarly by all functions. Now all functions
  are designed after McArdle & Anderson (*Ecology* 82, 290–297; 2001\).

* `dbrda` is a new function for distance-based Redundancy Analysis
  following McArdle & Anderson (*Ecology* 82, 290–297; 2001\). With
  metric dissimilarities, the function is equivalent to old `capscale`,
  but negative eigenvalues of semimetric indices are handled
  differently. In `dbrda` the dissimilarities are decomposed directly
  into conditions, constraints and residuals with their negative
  eigenvalues, and any of the components can have imaginary dimensions.
  Function is mostly compatible with `capscale` and other constrained
  ordination methods, but full compatibility cannot be achieved (see
  issue [\#140](https://github.com/vegandevs/vegan/issues/140) in
  Github). The function is based on the code by Pierre Legendre.

* The old `capscale` function for constrained ordination is still based
  only on real components, but the total inertia of the components is
  assessed similarly as in `dbrda`.
  
  The significance tests will differ from the previous version, but
  function `oldCapscale` will cast the `capscale` result to a similar
  form as previously.

* `adonis2` is a new function for permutational ANOVA of
  dissimilarities. It is based on the same algorithm as the `dbrda`. The
  function can perform overall tests of all independent variables as
  well as sequential and marginal tests of each term. The old `adonis`
  is still available, but it can only perform sequential tests. With
  same settings, `adonis` and `adonis2` give identical results (but see
  Github issue [\#156](https://github.com/vegandevs/vegan/issues/156)
  for differences).

* Function `varpart` can partition dissimilarities using the same
  algorithm as `dbrda`.

* Argument `sqrt.dist` takes square roots of dissimilarities and these
  can change many popular semimetric indices to metric distances in
  `capscale`, `dbrda`, `wcmdscale`, `adonis2`, `varpart` and
  `betadisper` (issue
  [\#179](https://github.com/vegandevs/vegan/issues/179) in Github).

* Lingoes and Cailliez adjustments change any dissimilarity into metric
  distance in `capscale`, `dbrda`, `adonis2`, `varpart`, `betadisper`
  and `wcmdscale`. Earlier we had only Cailliez adjustment in `capscale`
  (issue [\#179](https://github.com/vegandevs/vegan/issues/179) in
  Github).

* `RsquareAdj` works with `capscale` and `dbrda` and this allows using
  `ordiR2step` in model building.

## Bug Fixes

* `specaccum`: `plot` failed if line type (`lty`) was given. Reported by
  Lila Nath Sharma (Univ Bergen, Norway)

## New Functions

* `ordibar` is a new function to draw crosses of standard deviations or
  standard errors in ordination diagrams instead of corresponding
  ellipses.

* Several `permustats` results can be combined with a new `c()`
  function.

* New function `smbind` binds together null models by row, column or
  replication. If sequential models are bound together, they can be
  treated as parallel chains in subsequent analysis (e.g., after
  `as.mcmc`). See issue
  [\#164](https://github.com/vegandevs/vegan/issues/164) in Github.

## New Features

* Null model analysis was upgraded:
  
  New `"curveball"` algorithm provides a fast null model with fixed row
  and column sums for binary matrices after Strona et al. (*Nature
  Commun.* 5: 4114; 2014\).
  
  The `"quasiswap"` algorithm gained argument `thin` which can reduce
  the bias of null models.
  
  `"backtracking"` is now much faster, but it is still very slow, and
  provided mainly to allow comparison against better and faster methods.
  
  Compiled code can now be interrupted in null model simulations.

* `designdist` can now use beta diversity notation (`gamma`, `alpha`)
  for easier definition of beta diversity indices.

* `metaMDS` has new iteration strategy: Argument `try` gives the minimum
  number of random starts, and `trymax` the maximum number. Earlier we
  only hand `try` which gave the maximum number, but now we run at least
  `try` times. This reduces the risk of being trapped in a local optimum
  (issue [\#154](https://github.com/vegandevs/vegan/issues/154) in
  Github).
  
  If there were no convergent solutions, `metaMDS` will now tabulate
  stopping criteria (if `trace = TRUE`). This can help in deciding if
  any of the criteria should be made more stringent or the number of
  iterations increased. The documentation for `monoMDS` and `metaMDS`
  give more detailed information on convergence criteria.

* The `summary` of `permustats` prints now *P*-values, and the test
  direction (`alternative`) can be changed.
  
  The `qqmath` function of `permustats` can now plot standardized
  statistics. This is a partial solution to issue
  [\#172](https://github.com/vegandevs/vegan/issues/172) in Github.

* `MDSrotate` can rotate ordination to show maximum separation of factor
  levels (classes) using linear discriminant analysis (`lda` in
  [MASS](https://CRAN.R-project.org/package=MASS) package).

* `adipart`, `hiersimu` and `multipart` expose argument `method` to
  specify the null model.

* `RsquareAdj` works with `cca` and this allows using `ordiR2step` in
  model building. The code was developed by Dan McGlinn (issue
  [\#161](https://github.com/vegandevs/vegan/issues/161) in Github).
  However, `cca` still cannot be used in `varpart`.

* `ordiellipse` and `ordihull` allow setting colours, line types and
  other graphical parameters.
  
  The alpha channel can now be given also as a real number in 0 ... 1 in
  addition to integer 0 ... 255\.

* `ordiellipse` can now draw ellipsoid hulls that enclose points in a
  group.

* `ordicluster`, `ordisegments`, `ordispider` and `lines` and `plot`
  functions for `isomap` and `spantree` can use a mixture of colours of
  connected points. Their behaviour is similar as in analogous functions
  in the the [vegan3d](https://CRAN.R-project.org/package=vegan3d)
  package.

* `plot` of `betadisper` is more configurable. See issues
  [\#128](https://github.com/vegandevs/vegan/issues/128) and
  [\#166](https://github.com/vegandevs/vegan/issues/166) in Github for
  details.

* `text` and `points` methods for `orditkplot` respect stored graphical
  parameters.

* Environmental data for the Barro Colorado Island forest plots gained
  new variables from Harms et al. (*J. Ecol.* 89, 947–959; 2001\). Issue
  [\#178](https://github.com/vegandevs/vegan/issues/178) in Github.

## Deprecated and Defunct

* Function `metaMDSrotate` was removed and replaced with `MDSrotate`.

* `density` and `densityplot` methods for various vegan objects were
  deprecated and replaced with `density` and `densityplot` for
  `permustats`. Function `permustats` can extract the permutation and
  simulation results of vegan result objects.

# vegan 2.3-5

## Bug Fixes

* `eigenvals` fails with `prcomp` results in **R**-devel. The next
  version of `prcomp` will have an argument to limit the number of
  eigenvalues shown (`rank.`), and this breaks `eigenvals` in vegan.

* `calibrate` failed for `cca` and friends if `rank` was given.

# vegan 2.3-4

## Bug Fixes

* `betadiver` index `19` had wrong sign in one of its terms.

* `linestack` failed when the `labels` were given, but the input scores
  had no names. Reported by Jeff Wood (ANU, Canberra, ACT).

## Deprecated

* `vegandocs` is deprecated. Current **R** provides better tools for
  seeing extra documentation (`news()` and `browseVignettes()`).

## Vignettes

* All vignettes are built with standard **R** tools and can be browsed
  with `browseVignettes`. `FAQ-vegan` and `partitioning` were only
  accessible with `vegandocs` function.

## Building

* Dependence on external software `texi2dvi` was removed. Version 6\.1
  of `texi2dvi` was incompatible with **R** and prevented building
  vegan. The `FAQ-vegan` that was earlier built with `texi2dvi` uses now
  [knitr](https://CRAN.R-project.org/package=knitr). Because of this,
  vegan is now dependent on **R**-3\.0\.0\. Fixes issue
  [\#158](https://github.com/vegandevs/vegan/issues/158) in Github.

# vegan 2.3-3

## Bug Fixes

* `metaMDS` and `monoMDS` could fail if input dissimilarities were huge:
  in the reported case they were of magnitude 1E85\. Fixes issue
  [\#152](https://github.com/vegandevs/vegan/issues/152) in Github.

* Permutations failed if they were defined as
  [permute](https://CRAN.R-project.org/package=permute) control
  structures in `estaccum`, `ordiareatest`, `renyiaccum` and
  `tsallisaccum`. Reported by Dan Gafta (Cluj-Napoca) for `renyiaccum`.

* `rarefy` gave false warnings if input was a vector or a single
  sampling unit.

* Some extrapolated richness indices in `specpool` needed the number of
  doubletons (= number of species occurring in two sampling units), and
  these failed when only one sampling unit was supplied. The
  extrapolated richness cannot be estimated from a single sampling unit,
  but now such cases are handled smoothly instead of failing: observed
  non-extrapolated richness with zero standard error will be reported.
  The issue was reported in https://stackoverflow.com/questions/34027496/.

## New Features

* `treedist` and `treedive` refuse to handle trees with reversals, i.e,
  higher levels are more homogeneous than lower levels. Function
  `treeheight` will estimate their total height with absolute values of
  branch lengths. Function `treedive` refuses to handle trees with
  negative branch heights indicating negative dissimilarities. Function
  `treedive` is faster.

* `gdispweight` works when input data are in a matrix instead of a data
  frame.

* Input dissimilarities supplied in symmetric matrices or data frames
  are more robustly recognized by `anosim`, `bioenv` and `mrpp`.

# vegan 2.3-2

## Bug Fixes

* Printing details of a gridded permutation design would fail when the
  grid was at the within-plot level.

* `ordicluster` joined the branches at wrong coordinates in some cases.

* `ordiellipse` ignored weights when calculating standard errors (`kind
  = "se"`). This influenced plots of `cca`, and also influenced
  `ordiareatest`.

## New Features

* `adonis` and `capscale` functions recognize symmetric square matrices
  as dissimilarities. Formerly dissimilarities had to be given as
  `"dist"` objects such as produced by `dist` or `vegdist` functions,
  and data frames and matrices were regarded as observations x variables
  data which could confuse users (e.g., issue
  [\#147](https://github.com/vegandevs/vegan/issues/147)).

* `mso` accepts `"dist"` objects for the distances among locations as an
  alternative to coordinates of locations.

* `text`, `points` and `lines` functions for `procrustes` analysis
  gained new argument `truemean` which allows adding `procrustes` items
  to the plots of original analysis.

* `rrarefy` returns observed non-rarefied communities (with a warning)
  when users request subsamples that are larger than the observed
  community instead of failing. Function `drarefy` has been similar and
  returned sampling probabilities of 1, but now it also issues a
  warning. Fixes issue
  [\#144](https://github.com/vegandevs/vegan/issues/144) in Github.

# vegan 2.3-1

## Bug Fixes

* Permutation tests did not always correctly recognize ties with the
  observed statistic and this could result in too low `P`-values. This
  would happen in particular when all predictor variables were factors
  (classes). The changes concern functions `adonis`, `anosim`, `anova`
  and `permutest` functions for `cca`, `rda` and `capscale`, `permutest`
  for `betadisper`, `envfit`, `mantel` and `mantel.partial`, `mrpp`,
  `mso`, `oecosimu`, `ordiareatest`, `protest` and `simper`. This also
  fixes issues [\#120](https://github.com/vegandevs/vegan/issues/120)
  and [\#132](https://github.com/vegandevs/vegan/issues/132) in GitHub.

* Automated model building in constrained ordination (`cca`, `rda`,
  `capscale`) with `step`, `ordistep` and `ordiR2step` could fail if
  there were aliased candidate variables, or constraints that were
  completely explained by other variables already in the model. This was
  a regression introduced in vegan 2\.2-0\.

* Constrained ordination methods `cca`, `rda` and `capscale` treat
  character variables as factors in analysis, but did not return their
  centroids for plotting.

* Recovery of original data in `metaMDS` when computing WA scores for
  species would fail if the expression supplied to argument `comm` was
  long & got deparsed to multiple strings. `metaMDSdist` now returns the
  (possibly modified) data frame of community data `comm` as attribute
  `"comm"` of the returned `dist` object. `metaMDS` now uses this to
  compute the WA species scores for the NMDS. In addition, the deparsed
  expression for `comm` is now robust to long expressions. Reported by
  Richard Telford.

* `metaMDS` and `monoMDS` rejected dissimilarities with missing values.

* Function `rarecurve` did not check its input and this could cause
  confusing error messages. Now function checks that input data are
  integers that can be interpreted as counts on individuals and all
  sampling units have some species. Unchecked bad inputs were the reason
  for problems reported in https://stackoverflow.com/questions/30856909/.

## New Features and Functions

* Scaling of ordination axes in `cca`, `rda` and `capscale` can now be
  expressed with descriptive strings `"none"`, `"sites"`, `"species"` or
  `"symmetric"` to tell which kind of scores should be scaled by
  eigenvalues. These can be further modified with arguments `hill` in
  `cca` and `correlation` in `rda`. The old numeric scaling can still be
  used.

* The permutation data can be extracted from `anova` results of
  constrained ordination (`cca`, `rda`, `capscale`) and further analysed
  with `permustats` function.

* New data set `BCI.env` of site information for the Barro Colorado
  Island tree community data. Most useful variables are the UTM
  coordinates of sample plots. Other variables are constant or nearly
  constant and of little use in normal analysis.

# vegan 2.3-0

## Bug Fixes

* Constrained ordination functions `cca`, `rda` and `capscale` are now
  more robust. Scoping of data set names and variable names is much
  improved. This should fix numerous long-standing problems, for
  instance those reported by Benedicte Bachelot (in email) and Richard
  Telford (in Twitter), as well as issues
  [\#16](https://github.com/vegandevs/vegan/issues/16) and
  [\#100](https://github.com/vegandevs/vegan/issues/100) in GitHub.

* Ordination functions `cca` and `rda` silently accepted dissimilarities
  as input although their analysis makes no sense with these methods.
  Dissimilarities should be analysed with distance-based redundancy
  analysis (`capscale`).

* The variance of the conditional component was over-estimated in
  `goodness` of `rda` results, and results were wrong for partial RDA.
  The problems were reported in an
  [R-sig-ecology](https://stat.ethz.ch/pipermail/r-sig-ecology/2015-March/004936.html)
  message by Christoph von Redwitz.

## Windows

* `orditkplot` did not add file type identifier to saved graphics in
  Windows although that is required. The problem only concerned Windows
  OS.

## New Features and Functions

* `goodness` function for constrained ordination (`cca`, `rda`,
  `capscale`) was redesigned. Function gained argument `addprevious` to
  add the variation explained by previous ordination components to axes
  when `statistic = "explained"`. With this option, `model = "CCA"` will
  include the variation explained by partialled-out conditions, and
  `model = "CA"` will include the accumulated variation explained by
  conditions and constraints. The former behaviour was `addprevious =
  TRUE` for `model = "CCA"`, and `addprevious = FALSE` for `model =
  "CA"`. The argument will have no effect when `statistic = "distance"`,
  but this will always show the residual distance after all previous
  components. Formerly it displayed the residual distance only for the
  currently analysed model.

* Functions `ordiArrowMul` and `ordiArrowTextXY` are exported and can be
  used in normal interactive sessions. These functions are used to scale
  a bunch arrows to fit ordination graphics, and formerly they were
  internal functions used within other vegan functions.

* `orditkplot` can export graphics in SVG format. SVG is a vector
  graphics format which can be edited with several external programs,
  such as Illustrator and Inkscape.

* Rarefaction curve (`rarecurve`) and species accumulation models
  (`specaccum`, `fitspecaccum`) gained new functions to estimate the
  slope of curve at given location. Originally this was based on a
  response to an
  [R-SIG-ecology](https://stat.ethz.ch/pipermail/r-sig-ecology/2015-May/005038.html)
  query. For rarefaction curves, the function is `rareslope`, and for
  species accumulation models it is `specslope`.
  
  The functions are based on analytic equations, and can also be
  evaluated at interpolated non-integer values. In `specaccum` models
  the functions can be only evaluated for analytic models `"exact"`,
  `"rarefaction"` and `"coleman"`. With `"random"` and `"collector"`
  methods you can only use finite differences
  (`diff(fitted(<result.object>))`). Analytic functions for slope are
  used for all non-linear regression models known to `fitspecaccum`.

* Species accumulation models (`specaccum`) and non-liner regression
  models for species accumulation (`fitspecaccum`) work more
  consistently with weights. In all cases, the models are defined using
  the number of sites as independent variable, which with weights means
  that observations can be non-integer numbers of virtual sites. The
  `predict` models also use the number of sites with `newdata`, and for
  analytic models they can estimate the expected values for non-integer
  number of sites, and for non-analytic randomized or collector models
  they can interpolate on non-integer values.

* `fitspecaccum` gained support functions `AIC` and `deviance`.

* The `varpart` plots of four-component models were redesigned following
  Legendre, Borcard & Roberts *Ecology* 93, 1234–1240 (2012\), and they
  use now four ellipses instead of three circles and two rectangles. The
  components are now labelled in plots, and the circles and ellipses can
  be easily filled with transparent background colour.

# vegan 2.2-1

## General

* This is a maintenance release to avoid warning messages caused by
  changes in CRAN repository. The namespace usage is also more stringent
  to avoid warnings and notes in development versions of **R**.

## Installation

* vegan can be installed and loaded without tcltk package. The tcltk
  package is needed in `orditkplot` function for interactive editing of
  ordination graphics.

## Bug Fixes

* `ordisurf` failed if [gam](https://CRAN.R-project.org/package=gam)
  package was loaded due to namespace issues: some support functions of
  gam were used instead of
  [mgcv](https://CRAN.R-project.org/package=mgcv) functions.

* `tolerance` function failed for unconstrained correspondence analysis.

## New Features

* `estimateR` uses a more exact variance formula for bias-corrected Chao
  estimate of extrapolated number of species. The new formula may be
  unpublished, but it was derived following the guidelines of Chiu,
  Wang, Walther & Chao, *Biometrics* 70, 671–682 (2014\),
  [doi:10\.1111/biom.12200](https://doi.org/10.1111/biom.12200), online
  supplementary material.

* Diversity accumulation functions `specaccum`, `renyiaccum`,
  `tsallisaccum`, `poolaccum` and `estaccumR` use now
  [permute](https://CRAN.R-project.org/package=permute) package for
  permutations of the order of sampling sites. Normally these functions
  only need simple random permutation of sites, but restricted
  permutation of the permute package and user-supplied permutation
  matrices can be used.

* `estaccumR` function can use parallel processing.

* `linestack` accepts now expressions as labels. This allows using
  mathematical symbols and formula given as mathematical expressions.

# vegan 2.2-0

## General

* Several vegan functions can now use parallel processing for slow and
  repeating calculations. All these functions have argument `parallel`.
  The argument can be an integer giving the number of parallel
  processes. In unix-alikes (Mac OS, Linux) this will launch
  `"multicore"` processing and in Windows it will set up `"snow"`
  clusters as desribed in the documentation of the parallel package. If
  `option` `"mc.cores"` is set to an integer \> 1, this will be used to
  automatically start parallel processing. Finally, the argument can
  also be a previously set up `"snow"` cluster which will be used both
  in Windows and in unix-alikes. Vegan vignette on Design decision
  explains the implementation (use `vegandocs("decission")`, and
  parallel package has more extensive documentation on parallel
  processing in **R**.
  
  The following function use parallel processing in analysing
  permutation statistics: `adonis`, `anosim`, `anova.cca` (and
  `permutest.cca`), `mantel` (and `mantel.partial`), `mrpp`,
  `ordiareatest`, `permutest.betadisper` and `simper`. In addition,
  `bioenv` can compare several candidate sets of models in paralle,
  `metaMDS` can launch several random starts in parallel, and `oecosimu`
  can evaluate test statistics for several null models in parallel.

* All permutation tests are based on the
  [permute](https://CRAN.R-project.org/package=permute) package which
  offers strong tools for restricted permutation. All these functions
  have argument `permutations`. The default usage of simple
  non-restricted permutations is achieved by giving a single integer
  number. Restricted permutations can be defined using the `how`
  function of the permute package. Finally, the argument can be a
  permutation matrix where rows define permutations. It is possible to
  use external or user constructed permutations.
  
  See `help(permutations)` for a brief introduction on permutations in
  vegan, and permute package for the full documention. The vignette of
  the permute package can be read from vegan with command
  `vegandocs("permutations")`.
  
  The following functions use the
  [permute](https://CRAN.R-project.org/package=permute) package:
  `CCorA`, `adonis`, `anosim`, `anova.cca` (plus associated
  `permutest.cca`, `add1.cca`, `drop1.cca`, `ordistep`, `ordiR2step`),
  `envfit` (plus associated `factorfit` and `vectorfit`), `mantel` (and
  `mantel.partial`), `mrpp`, `mso`, `ordiareatest`,
  `permutest.betadisper`, `protest` and `simper`.

* Community null model generation has been completely redesigned and
  rewritten. The communities are constructed with new `nullmodel`
  function and defined in a low level `commsim` function. The actual
  null models are generated with a `simulate` function that builds an
  array of null models. The new null models include a wide array of
  quantitative models in addition to the old binary models, and users
  can plug in their own generating functions. The basic tool invoking
  and analysing null models is `oecosimu`. The null models are often
  used only for the analysis of nestedness, but the implementation in
  `oecosimu` allows analysing any statistic, and null models are better
  seen as an alternative to permutation tests.

## Installation

* vegan package dependencies and namespace imports were adapted to
  changes in **R**, and no more trigger warnings and notes in package
  tests.

* Three-dimensional ordination graphics using
  [scatterplot3d](https://CRAN.R-project.org/package=scatterplot3d) for
  static plots and [rgl](https://CRAN.R-project.org/package=rgl) for
  dynamic plots were removed from vegan and moved to a companion package
  [vegan3d](https://CRAN.R-project.org/package=vegan3d). The package is
  available in CRAN.

## New Functions

* Function `dispweight` implements dispersion weighting of Clarke et al.
  (*Marine Ecology Progress Series*, 320, 11–27\). In addition, we
  implemented a new method for generalized dispersion weighting
  `gdispweight`. Both methods downweight species that are significantly
  over-dispersed.

* New `hclust` support functions `reorder`, `rev` and `scores`.
  Functions `reorder` and `rev` are similar as these functions for
  `dendrogram` objects in base **R**. However, `reorder` can use (and
  defaults to) weighted mean. In weighted mean the node average is
  always the mean of member leaves, whereas the `dendrogram` uses always
  unweighted means of joined branches.

* Function `ordiareatest` supplements `ordihull` and `ordiellipse` and
  provides a randomization test for the one-sided alternative hypothesis
  that convex hulls or ellipses in two-dimensional ordination space have
  smaller areas than with randomized groups.

* Function `permustats` extracts and inspects permutation results with
  support functions `summary`, `density`, `densityplot`, `qqnorm` and
  `qqmath`. The `density` and `qqnorm` are standard **R** tools that
  only work with one statistic, and `densityplot` and `qqmath` are
  lattice graphics that work with univariate and multivariate
  statistics. The results of following functions can be extracted:
  `anosim`, `adonis`, `mantel` (and `mantel.partial`), `mrpp`,
  `oecosimu`, `permustest.cca` (but not the corresponding `anova`
  methods), `permutest.betadisper`, and `protest`.

* `stressplot` functions display the ordination distances at given
  number of dimensions against original distances. The method functins
  are similar to `stressplot` for `metaMDS`, and always use the inherent
  distances of each ordination method. The functions are available for
  the results `capscale`, `cca`, `princomp`, `prcomp`, `rda`, and
  `wcmdscale`.

## Bug Fixes

* `cascadeKM` of only one group will be `NA` instead of a random value.

* `ordiellipse` can handle points exactly on a line, including only two
  points (with a warning).

* plotting `radfit` results for several species failed if any of the
  communities had no species or had only one species.

* `RsquareAdj` for `capscale` with negative eigenvalues will now report
  `NA` instead of using biased method of `rda` results.

* `simper` failed when a group had only a single member.

## New Features

* `anova.cca` functions were re-written to use the permute package. Old
  results may not be exactly reproduced, and models with missing data
  may fail in several cases. There is a new option of analysing a
  sequence of models against each other.

* `simulate` functions for `cca` and `rda` can return several
  simulations in a `nullmodel` compatible object. The functions can
  produce simulations with correlated errors (also for `capscale`) in
  parametric simulation with Gaussian error.

* `bioenv` can use Manhattan, Gower and Mahalanobis distances in
  addition to the default Euclidean. New helper function `bioenvdist`
  can extract the dissimilarities applied in best model or any other
  model.

* `metaMDS(..., trace = 2)` will show convergence information with the
  default `monoMDS` engine.

* Function `MDSrotate` can rotate a `k`-dimensional ordination to `k-1`
  variables. When these variables are correlated (like usually is the
  case), the vectors can also be correlated to previously rotated
  dimensions, but will be uncorrelated to all later ones.

* vegan 2\.0-10 changed the weighted `nestednodf` so that weighted
  analysis of binary data was equivalent to binary analysis. However,
  this broke the equivalence to the original method. Now the function
  has an argument `wbinary` to select the method of analysis. The
  problem was reported and a fix submitted by Vanderlei Debastiani
  (Universidade Federal do Rio Grande do Sul, Brasil).

* `ordiellipse`, `ordihull` and `ordiellipse` can handle missing values
  in `groups`.

* `ordispider` can now use spatial medians instead of means.

* `rankindex` can use Manhattan, Gower and Mahalanobis distance in
  addition to the default Euclidean.

* User can set colours and line types in function `rarecurve` for
  plotting rarefaction curves.

* `spantree` gained a support function `as.hclust` to change the minimum
  spanning tree into an `hclust` tree.

* `fitspecaccum` can do weighted analysis. Gained `lines` method.

* Functions for extrapolated number of species or for the size of
  species pool using Chao method were modified following Chiu et al.,
  *Biometrics* 70, 671–682 (2014\).
  
  Incidence based `specpool` can now use (and defaults to) small sample
  correction with number of sites as the sample size. Function uses
  basic Chao extrapolation based on the ratio of singletons and
  doubletons, but switches now to bias corrected Chao extrapolation if
  there are no doubletons (species found twice). The variance formula
  for bias corrected Chao was derived following the supporting on line
  material of
  [doi:10\.1111/biom.12200](https://doi.org/10.1111/biom.12200) and
  differs slightly from Chiu et al. (2014\).
  
  The `poolaccum` function was changed similarly, but the small sample
  correction is used always.
  
  The abundance based `estimateR` uses bias corrected Chao
  extrapolation, but earlier it estimated its variance with classic Chao
  model. Now we use the widespread approximate estimate from EstimateS
  for variance.
  
  With these changes these functions are more similar to **EstimateS**

* `tabasco` uses now `reorder.hclust` for `hclust` object for better
  ordering than previously when it cast trees to `dendrogram` objects.

* `treedive` and `treedist` default now to `match.force = TRUE` and can
  be silenced with `verbose = FALSE`.

* `vegdist` gained Mahalanobis distance.

* Nomenclature updated in plant community data with the help of
  Taxonstand and taxize packages. The taxonomy of the `dune` data was
  adapted to the same sources and APG III. `varespec` and `dune` use
  8-character names (4 from genus + 4 from species epithet). New data
  set on phylogenetic distances for `dune` was extracted from Zanne et
  al. (*Nature* 506, 89–92; 2014\).

* User configurable plots for `rarecurve`.

## Deprecated and Defunct

* `strata` are deprecated in permutations. It is still accepted but will
  be phased out in next releases. Use `how` of permute package.

* `cca`, `rda` and `capscale` do not return scores scaled by
  eigenvalues: use `scores` function to extract scaled results.

* `commsimulator` is deprecated. Replace `commsimulator(x, method)` with
  `simulate(nullmodel(x, method))`.

* `density` and `densityplot` for permutation results are deprecated:
  use `permustats` with its `density` and `densityplot` method.

# vegan 2.0-10

## General

* This version is adapted to the changes in permute package version
  0\.8-0 and no more triggers NOTEs in package checks. This release may
  be the last of the 2\.0 series, and the next vegan release is
  scheduled to be a major release with newly designed `oecosimu` and
  community pattern simulation, support for parallel processing, and
  full support of the permute package. If you are interested in these
  developments, you may try the development versions of vegan in
  [GitHub](https://github.com/jarioksa/vegan) and report the problems
  and user experience to us.

## Bug Fixes

* `envfit` function assumed that all external variables were either
  numeric or factors, and failed if they were, say, character strings.
  Now only numeric variables are taken as continuous vectors, and all
  other variables (character strings, logical) are coerced to factors if
  possible. The function also should work with degenerate data, like
  only one level of a factor or a constant value of a continuous
  environmental variable. The ties were wrongly in assessing permutation
  `P`-values in `vectorfit`.

* `nestednodf` with quantitative data was not consistent with binary
  models, and the fill was wrongly calculated with quantitative data.

* `oecosimu` now correctly adapts displayed quantiles of simulated
  values to the `alternative` test direction.

* `renyiaccum` plotting failed if only one level of diversity `scale`
  was used.

## New Features

* The Kempton and Taylor algorithm was found unreliable in `fisherfit`
  and `fisher.alpha`, and now the estimation of Fisher &alpha; is only
  based on the number of species and the number of individuals. The
  estimation of standard errors and profile confidence intervals also
  had to be scrapped.

* `renyiaccum`, `specaccum` and `tsallisaccum` functions gained `subset`
  argument.

* `renyiaccum` can now add a `collector` curve to to the analysis. The
  collector curve is the diversity accumulation in the order of the
  sampling units. With an interesting ordering or sampling units this
  allows comparing actual species accumulations with the expected
  randomized accumulation.

* `specaccum` can now perform weighted accumulation using the sampling
  effort as weights.

# vegan 2.0-9

* This version is released due to changes in programming interface and
  testing procedures in **R** 3\.0\.2\. If you are using an older
  version of **R**, there is no need to upgrade vegan. There are no new
  features nor bug fixes. The only user-visible changes are in
  documentation and in output messages and formatting. Because of **R**
  changes, this version is dependent on **R** version 2\.14\.0 or newer
  and on lattice package.

# vegan 2.0-8

## General

* This is a maintenance release that fixes some issues raised by changed
  in **R** toolset for processing vignettes. In the same we also fix
  some typographic issues in the vignettes.

## New Features

* `ordisurf` gained new arguments for more flexible definition of fitted
  models to better utilize the mgcv`::gam` function.
  
  The linewidth of contours can now be set with the argument `lwd`.

* Labels to arrows are positioned in a better way in `plot` functions
  for the results of `envfit`, `cca`, `rda` and `capscale`. The labels
  should no longer overlap the arrow tips.

* The setting test direction is clearer in `oecosimu`.

* `ordipointlabel` gained a `plot` method that can be used to replot the
  saved result.

# vegan 2.0-7

## New Functions

* `tabasco()` is a new function for graphical display of community data
  matrix. Technically it is an interface to **R** `heatmap`, but its use
  is closer to vegan function `vegemite`. The function can reorder the
  community data matrix similarly as `vegemite`, for instance, by
  ordination results. Unlike `heatmap`, it only displays dendrograms if
  supplied by the user, and it defaults to re-order the dendrograms by
  correspondence analysis. Species are ordered to match site ordering or
  like determined by the user.

## Bug Fixes

* Function `fitspecaccum(..., model = "asymp")` fitted logistic model
  instead of asymptotic model (or the same as `model = "logis"`).

* `nestedtemp()` failed with very sparse data (fill `< 0.38`%).

## New Features

* The `plot` function for constrained ordination results (`cca`, `rda`,
  `capscale`) gained argument `axis.bp` (defaults `TRUE`) which can be
  used to suppress axis scale for biplot arrays.

* Number of iterations in nonmetric multidimensional scaling (NMDS) can
  be set with keyword `maxit` (defaults `200`) in `metaMDS`.

## Deprecated

* The result objects of `cca`, `rda` and `capscale` will no longer have
  scores `u.eig`, `v.eig` and `wa.eig` in the future versions of vegan.
  This change does not influence normal usage, because vegan functions
  do not need these items. However, external scripts and packages may
  need changes in the future versions of vegan.

# vegan 2.0-6

## Bug Fixes

* The species scores were scaled wrongly in `capscale()`. They were
  scaled correctly only when Euclidean distances were used, but usually
  `capscale()` is used with non-Euclidean distances. Most graphics will
  change and should be redone. The change of scaling mainly influences
  the spread of species scores with respect to the site scores.

* Function `clamtest()` failed to set the minimum abundance threshold in
  some cases. In addition, the output was wrong when some of the
  possible species groups were missing. Both problems were reported by
  Richard Telford (Bergen, Norway).

* Plotting an object fitted by `envfit()` would fail if `p.max` was used
  and there were unused levels for one or more factors. The unused
  levels could result from deletion of observations with missing values
  or simply as the result of supplying a subset of a larger data set to
  `envfit()`.

* `multipart()` printed wrong information about the analysis type (but
  did the analysis correctly). Reported by Valerie Coudrain.

* `oecosimu()` failed if its `nestedfun` returned a data frame. A more
  fundamental fix will be in vegan 2\.2-0, where the structure of the
  `oecosimu()` result will change.

* The plot of two-dimensional `procrustes()` solutions often draw
  original axes in a wrong angle. The problem was reported by Elizabeth
  Ottesen (MIT).

* Function `treedive()` for functional or phylogenetic diversity did not
  correctly match the species names between the community data and
  species tree when the tree contained species that did not occur in the
  data. Related function `treedist()` for phylogenetic distances did not
  try to match the names at all.

## New Features

* The output of `capscale()` displays the value of the additive constant
  when argument `add = TRUE` was used.

* `fitted()` functions for `cca()`, `rda()` and `capscale()` can now
  return conditioned (partial) component of the response: Argument
  `model` gained a new alternative `model = "pCCA"`.

* `dispindmorisita()` output gained a new column for Chi-squared based
  probabilities that the null hypothesis (random distribution) is true.

* `metaMDS()` and `monoMDS()` have new default convergence criteria.
  Most importantly, scale factor of the gradient (`sfgrmin`) is
  stricter. The former limit was too slack with large data sets and
  iterations stopped early without getting close to the solution. In
  addition, `scores()` ignore now requests to dimensions beyond those
  calculated instead of failing, and `scores()` for `metaMDS()` results
  do not drop dimensions.

* `msoplot()` gained `legend` argument for positioning the legend.

* Nestedness function `nestednodf()` gained a `plot` method.

* `ordiR2step()` gained new argument `R2scope` (defaults `TRUE`) which
  can be used to turn off the criterion of stopping when the adjusted
  _R_<sup>2</sup> of the current model exceeds that of the scope. This option
  allows model building when the `scope` would be overdetermined (number
  of predictors higher than number of observations).
  
  `ordiR2step()` now handles partial redundancy analysis (pRDA).

* `orditorp()` gained argument `select` to select the rows or columns of
  the results to display.

* `protest()` prints the standardized residual statistic squared m12
  in addition to the squared Procrustes correlation _R_<sup>2</sup>. Both
  were calculated, but only the latter was displayed.
  
  Permutation tests are much faster in `protest()`. Instead of calling
  repeatedly `procrustes()`, the goodness of fit statistic is evaluated
  within the function.

* `wcmdscale()` gained methods for `print`, `plot` etc. of the results.
  These methods are only used if the full `wcmdscale` result is returned
  with, e.g., argument `eig = TRUE`. The default is still to return only
  a matrix of scores similarly as the standard **R** function
  `cmdscale()`, and in that case the new methods are not used.

# vegan 2.0-5

## Bug Fixes

* `anova(<cca_object>, ...)` failed with `by = "axis"` and `by =
  "term"`. The bug was reported by Dr Sven Neulinger (Christian Albrecht
  University, Kiel, Germany).

* `radlattice` did not honour argument `BIC = TRUE`, but always
  displayed AIC.

## New Functions

* Most vegan functions with permutation tests have now a `density`
  method that can be used to find empirical probability distributions of
  permutations. There is a new `plot` method for these functions that
  displays both the density and the observed statistic. The `density`
  function is available for `adonis`, `anosim`, `mantel`,
  `mantel.partial`, `mrpp`, `permutest.cca` and `procrustes`.
  
  Function `adonis` can return several statistics, and it has now a
  `densityplot` method (based on lattice).
  
  Function `oecosimu` already had `density` and `densityplot`, but they
  are now similar to other vegan methods, and also work with `adipart`,
  `hiersimu` and `multipart`.

* `radfit` functions got a `predict` method that also accepts arguments
  `newdata` and `total` for new ranks and site totals for prediction.
  The functions can also interpolate to non-integer “ranks”, and in some
  models also extrapolate.

## New Features

* Labels can now be set in the `plot` of `envfit` results. The labels
  must be given in the same order that the function uses internally, and
  new support function `labels` can be used to display the default
  labels in their correct order.

* Mantel tests (functions `mantel` and `mantel.partial`) gained argument
  `na.rm` which can be used to remove missing values. This options
  should be used with care: Permutation tests can be biased if the
  missing values were originally in matching or fixed positions.

* `radfit` results can be consistently accessed with the same methods
  whether they were a single model for a single site, all models for a
  single site or all models for all sites in the data. All functions now
  have methods `AIC`, `coef`, `deviance`, `logLik`, `fitted`, `predict`
  and `residuals`.

## Installation and Building

* Building of vegan vignettes failed with the latest version of LaTeX
  (TeXLive 2012\).

* **R** versions later than 2\.15-1 (including development version)
  report warnings and errors when installing and checking vegan, and you
  must upgrade vegan to this version. The warnings concern functions
  `cIndexKM` and `betadisper`, and the error occurs in `betadisper`.
  These errors and warnings were triggered by internal changes in **R**.

# vegan 2.0-4

## Bug Fixes

* `adipart` assumed constant gamma diversity in simulations when
  assessing the `P`-value. This could give biased results if the null
  model produces variable gamma diversities and option `weights =
  "prop"` is used. The default null model (`"r2dtable"`) and the default
  option (`weights = "unif"`) were analysed correctly.

* `anova(<prc-object>, by = "axis")` and other `by` cases failed due to
  ‘NAMESPACE’ issues.

* `clamtest` wrongly used frequencies instead of the counts when
  calculating sample coverage. No detectable differences were produced
  when rerunning examples from Chazdon et al. 2011 and vegan help page.

* `envfit` failed with unused factor levels.

* `predict` for `cca` results with `type = "response"` or `type =
  "working"` failed with `newdata` if the number of rows did not match
  with the original data. Now the `newdata` is ignored if it has a wrong
  number of rows. The number of rows must match because the results in
  `cca` must be weighted by original row totals. The problem did not
  concern `rda` or `capscale` results which do not need row weights.
  Reported by Glenn De'ath.

## New Features

* Functions for diversity partitioning (`adipart`, `hiersimu` and
  `multipart`) have now `formula` and `default` methods. The `formula`
  method is identical to the previous functions, but the `default`
  method can take two matrices as input.
  
  Functions `adipart` and `multipart` can be used for fast and easy
  overall partitioning to alpha, beta and gamma diversities by omitting
  the argument describing the hierarchy.

* The method in `betadisper` is biased with small sample sizes. The
  effects of the bias are strongest with unequal sample sizes. A bias
  adjusted version was developed by Adrian Stier and Ben Bolker, and can
  be invoked with argument `bias.adjust` (defaults to `FALSE`).

* `bioenv` accepts dissimilarities (or square matrices that can be
  interpreted as dissimilarities) as an alternative to community data.
  This allows using other dissimilarities than those available in
  `vegdist`.

* `plot` function for `envfit` results gained new argument `bg` that can
  be used to set background colour for plotted labels.

* `msoplot` is more configurable, and allows, for instance, setting
  y-axis limits.

* Hulls and ellipses are now filled using semitransparent colours in
  `ordihull` and `ordiellipse`, and the user can set the degree of
  transparency with a new argument `alpha`. The filled shapes are used
  when these functions are called with argument `draw = "polygon"`.
  Function `ordihull` puts labels (with argument `label = TRUE`) now in
  the real polygon centre.

* `ordiplot3d` returns function `envfit.convert` and the projected
  location of the `origin`. Together these can be used to add `envfit`
  results to existing `ordiplot3d` plots.
  
  Equal aspect ratio cannot be set exactly in `ordiplot3d` because
  underlying core routines do not allow this. Now `ordiplot3d` sets
  equal axis ranges, and the documents urge users to verify that the
  aspect ratio is reasonably equal and the graph looks like a cube. If
  the problems cannot be solved in the future, `ordiplot3d` may be
  removed from next releases of vegan.

* Function `ordipointlabel` gained argument to `select` only some of the
  items for plotting. The argument can be used only with one set of
  points.

# vegan 2.0-3

## New Functions

* Added new nestedness functions `nestedbetasor` and `nestedbetajac`
  that implement multiple-site dissimilarity indices and their
  decomposition into turnover and nestedness components following
  Baselga (*Global Ecology and Biogeography* 19, 134–143; 2010\).

* Added function `rarecurve` to draw rarefaction curves for each row
  (sampling unit) of the input data, optionally with lines showing
  rarefied species richness with given sample size for each curve.

* Added function `simper` that implements “similarity percentages” of
  Clarke (*Australian Journal of Ecology* 18, 117–143; 1993\). The
  method compares two or more groups and decomposes the average
  between-group Bray-Curtis dissimilarity index to contributions by
  individual species. The code was developed in
  [GitHub](https://github.com/jarioksa/vegan) by Eduard Szöcs (Uni
  Landau, Germany).

## Bug Fixes

* `betadisper()` failed when the `groups` was a factor with empty
  levels.

* Some constrained ordination methods and their support functions are
  more robust in border cases (completely aliased effects, saturated
  models, user requests for non-existng scores etc). Concerns
  `capscale`, `ordistep`, `varpart`, `plot` function for constrained
  ordination, and `anova(<cca.object>, by = "margin")`.

* The `scores` function for `monoMDS` did not honour `choices` argument
  and hence dimensions could not be chosen in `plot`.

* The default `scores` method failed if the number of requested axes was
  higher than the ordination object had. This was reported as an error
  in `ordiplot` in
  [R-sig-ecology](https://stat.ethz.ch/pipermail/r-sig-ecology/2012-February/002768.html)
  mailing list.

## New Features

* `metaMDS` argument `noshare = 0` is now regarded as a numeric
  threshold that always triggers extended dissimilarities
  (`stepacross`), instead of being treated as synonymous with `noshare =
  FALSE` which always suppresses extended dissimilarities.

* Nestedness discrepancy index `nesteddisc` gained a new argument that
  allows user to set the number of iterations in optimizing the index.

* `oecosimu` displays the mean of simulations and describes alternative
  hypothesis more clearly in the printed output.

* Implemented adjusted _R_<sup>2</sup> for partial RDA. For partial model `rda(Y ~
  X1 + Condition(X2))` this is the same as the component `[a] = X1|X2`
  in variance partition in `varpart` and describes the marginal (unique)
  effect of constraining term to adjusted _R_<sup>2</sup>.

* Added Cao dissimilarity (CYd) as a new dissimilarity method in
  `vegdist` following Cao et al., *Water Envir Res* 69, 95–106 (1997\).
  The index should be good for data with high beta diversity and
  variable sampling intensity. Thanks to consultation to Yong Cao (Univ
  Illinois, USA).

# vegan 2.0-2

## Bug Fixes

* Function `capscale` failed if constrained component had zero rank.
  This happened most likely in partial models when the conditions
  aliased constraints. The problem was observed in `anova(..., by
  ="margin")` which uses partial models to analyses the marginal
  effects, and was reported in an email message to [R-News mailing
  list](https://stat.ethz.ch/pipermail/r-help/2011-October/293077.html).

* `stressplot` and `goodness` sometimes failed when `metaMDS` was based
  on `isoMDS` (MASS package) because `metaMDSdist` did not use the same
  defaults for step-across (extended) dissimilarities as `metaMDS(...,
  engine = "isoMDS")`. The change of defaults can also influence
  triggering of step-across in `capscale(..., metaMDSdist = TRUE)`.

* `adonis` contained a minor bug resulting from incomplete
  implementation of a speed-up that did not affect the results. In
  fixing this bug, a further bug was identified in transposing the hat
  matrices. This second bug was only active following fixing of the
  first bug. In fixing both bugs, a speed-up in the internal f.test()
  function is fully realised. Reported by Nicholas Lewin-Koh.

## New Features

* `ordiarrows` and `ordisegments` gained argument `order.by` that gives
  a variable to sort points within `groups`. Earlier the points were
  assumed to be in order.

* Function `ordispider` invisibly returns the coordinates to which the
  points were connected. Typically these are class centroids of each
  point, but for constrained ordination with no `groups` they are the LC
  scores.

# vegan 2.0-1

## New Functions

* `clamtest`: new function to classify species as generalists and
  specialists in two distinct habitats (CLAM test of Chazdon et al.,
  *Ecology* 92, 1332–1343; 2011\). The test is based on multinomial
  distribution of individuals in two habitat types or sampling units,
  and it is applicable only to count data with no over-dispersion.

* `as.preston` gained `plot` and `lines` methods, and `as.fisher` gained
  `plot` method (which also can add items to existing plots). These are
  similar as `plot` and `lines` for `prestonfit` and `fisherfit`, but
  display only data without the fitted lines.

* `raupcrick`: new function to implement Raup-Crick dissimilarity as a
  probability of number of co-occurring species with occurrence
  probabilities proportional to species frequencies. Vegan has
  Raup-Crick index as a choice in `vegdist`, but that uses equal
  sampling probabilities for species and analytic equations. The new
  `raupcrick` function uses simulation with `oecosimu`. The function
  follows Chase et al. (2011\) *Ecosphere* 2:art24
  \[[doi:10\.1890/ES10-00117\.1](https://doi.org/10.1890/ES10-00117.1)\],
  and was developed with the consultation of Brian Inouye.

## Bug Fixes

* Function `meandist` could scramble items and give wrong results,
  especially when the `grouping` was numerical. The problem was reported
  by Dr Miguel Alvarez (Univ. Bonn).

* `metaMDS` did not reset `tries` when a new model was started with a
  `previous.best` solution from a different model.

* Function `permatswap` for community null models using quantitative
  swap never swapped items in a 2x2 submatrix if all cells were
  filled.

* The result from `permutest.cca` could not be `update`d because of a
  ‘NAMESPACE’ issue.

* **R** 2\.14\.0 changed so that it does not accept using `sd()`
  function for matrices (which was the behaviour at least since **R**
  1\.0-0\), and several vegan functions were changed to adapt to this
  change (`rda`, `capscale`, `simulate` methods for `rda`, `cca` and
  `capscale`). The change in **R** 2\.14\.0 does not influence the
  results but you probably wish to upgrade vegan to avoid annoying
  warnings.

## Analyses

* `nesteddisc` is slacker and hence faster when trying to optimize the
  statistic for tied column frequencies. Tracing showed that in most
  cases an improved ordering was found rather early in tries, and the
  results are equally good in most cases.

# vegan 2.0-0

## General

* Peter Minchin joins the vegan team.

* vegan implements standard **R** ‘NAMESPACE’. In general, `S3` methods
  are not exported which means that you cannot directly use or see
  contents of functions like `cca.default`, `plot.cca` or
  `anova.ccabyterm`. To use these functions you should rely on **R**
  delegation and simply use `cca` and for its result objects use `plot`
  and `anova` without suffix `.cca`. To see the contents of the function
  you can use `:::`, such as `vegan:::cca.default`. This change may
  break packages, documents or scripts that rely on non-exported names.

* vegan depends on the permute package. This package provides powerful
  tools for restricted permutation schemes. All vegan permutation will
  gradually move to use permute, but currently only `betadisper` uses
  the new feature.

## New Functions

* `monoMDS`: a new function for non-metric multidimensional scaling
  (NMDS). This function replaces `MASS::isoMDS` as the default method in
  `metaMDS`. Major advantages of `monoMDS` are that it has ‘weak’
  (‘primary’) tie treatment which means that it can split tied
  observed dissimilarities. ‘Weak’ tie treatment improves ordination of
  heterogeneous data sets, because maximum dissimilarities of `1` can be
  split. In addition to global NMDS, `monoMDS` can perform local and
  hybrid NMDS and metric MDS. It can also handle missing and zero
  dissimilarities. Moreover, `monoMDS` is faster than previous
  alternatives. The function uses `Fortran` code written by Peter
  Minchin.

* `MDSrotate` a new function to replace `metaMDSrotate`. This function
  can rotate both `metaMDS` and `monoMDS` results so that the first axis
  is parallel to an environmental vector.

* `eventstar` finds the minimum of the evenness profile on the Tsallis
  entropy, and uses this to find the corresponding values of diversity,
  evenness and numbers equivalent following Mendes et al. (*Ecography*
  31, 450-456; 2008\). The code was contributed by Eduardo Ribeira Cunha
  and Heloisa Beatriz Antoniazi Evangelista and adapted to vegan by
  Peter Solymos.

* `fitspecaccum` fits non-linear regression models to the species
  accumulation results from `specaccum`. The function can use new
  self-starting species accumulation models in vegan or other
  self-starting non-linear regression models in **R**. The function can
  fit Arrhenius, Gleason, Gitay, Lomolino (in vegan), asymptotic,
  Gompertz, Michaelis-Menten, logistic and Weibull (in base **R**)
  models. The function has `plot` and `predict` methods.

* Self-starting non-linear species accumulation models `SSarrhenius`,
  `SSgleason`, `SSgitay` and `SSlomolino`. These can be used with
  `fitspecaccum` or directly in non-linear regression with `nls`. These
  functions were implemented because they were found good for
  species-area models by Dengler (*J. Biogeogr.* 36, 728-744; 2009\).

## New Features

* `adonis`, `anosim`, `meandist` and `mrpp` warn on negative
  dissimilarities, and `betadisper` refuses to analyse them. All these
  functions expect dissimilarities, and giving something else (like
  correlations) probably is a user error.

* `betadisper` uses restricted permutation of the permute package.

* `metaMDS` uses `monoMDS` as its default ordination engine. Function
  gains new argument `engine` that can be used to alternatively select
  `MASS::isoMDS`. The default is not to use `stepacross` with `monoMDS`
  because its ‘weak’ tie treatment can cope with tied maximum
  dissimilarities of one. However, `stepacross` is the default with
  `isoMDS` because it cannot handle adequately these tied maximum
  dissimilarities.

* `specaccum` gained `predict` method which uses either linear or spline
  interpolation for data between observed points. Extrapolation is
  possible with spline interpolation, but may make little sense.

* `specpool` can handle missing values or empty factor levels in the
  grouping factor `pool`. Now also checks that the length of the `pool`
  matches the number of observations.

## Deprecated and Defunct

* `metaMDSrotate` was replaced with `MDSrotate` that can also handle the
  results of `monoMDS`.

* `permuted.index2` and other “new” permutation code was removed in
  favour of the permute package. This code was not intended for normal
  use, but packages depending on that code in vegan should instead
  depend on permute.

## Analyses

* `treeheight` uses much snappier code. The results should be unchanged.
