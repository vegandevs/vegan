\name{envfit}
\alias{envfit}
\alias{envfit.default}
\alias{envfit.formula}
\alias{vectorfit}
\alias{factorfit}
\alias{plot.envfit}
\alias{scores.envfit}
\alias{labels.envfit}

\title{Fits an Environmental Vector or Factor onto an Ordination }

\description{

  The function fits environmental vectors or factors onto an
  ordination. The projections of points onto vectors have maximum
  correlation with corresponding environmental variables, and the
  factors show the averages of factor levels. For continuous varaibles
  this is equal to fitting a linear trend surface (plane in 2D) for a
  variable (see \code{\link{ordisurf}}); this trend surface can be
  presented by showing its gradient (direction of steepest increase)
  using an arrow. The environmental variables are the dependent
  variables that are explained by the ordination scores, and each
  dependent variable is analysed separately.

}

\usage{
\method{envfit}{default}(ord, env, permutations = 999, strata = NULL, 
   choices=c(1,2),  display = "sites", w, na.rm = FALSE, ...)
\method{envfit}{formula}(formula, data, ...)
\method{plot}{envfit}(x, choices = c(1,2), labels, arrow.mul, at = c(0,0), 
   axis = FALSE, p.max = NULL, col = "blue", bg, add = TRUE, ...)
\method{scores}{envfit}(x, display, choices, arrow.mul=1, tidy = FALSE, ...)
vectorfit(X, P, permutations = 0, strata = NULL, w, ...)
factorfit(X, P, permutations = 0, strata = NULL, w, ...)
}

\arguments{
  \item{ord}{An ordination object or other structure from which the
    ordination \code{\link{scores}} can be extracted (including a data
    frame or matrix of scores).}
  \item{env}{Data frame, matrix or vector of environmental
    variables. The variables can be of mixed type (factors, continuous
    variables) in data frames.}
  \item{X}{Matrix or data frame of ordination scores.}
  \item{P}{Data frame, matrix or vector of environmental
    variable(s). These must be continuous for \code{vectorfit} and
    factors or characters for \code{factorfit}. }
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices. Set \code{permutations = 0} to skip
    permutations.}
  \item{formula, data}{Model  \code{\link{formula}} and data.  }
  \item{na.rm}{Remove points with missing values in ordination scores
    or environmental variables. The operation is casewise: the whole
    row of data is removed if there is a missing value and 
    \code{na.rm = TRUE}.}
    
  \item{x}{A result object from \code{envfit}. For \code{ordiArrowMul}
    and \code{ordiArrowTextXY} this must be a two-column matrix (or
    matrix-like object) containing the coordinates of arrow heads on
    the two plot axes, and other methods extract such a structure from
    the \code{envfit} results.}
    
  \item{choices}{Axes to plotted.}
  \item{tidy}{Return scores that are compatible with \CRANpkg{ggplot2}:
    all scores are in a single \code{data.frame}, score type is
    identified by factor variable \code{scores} (\code{"vectors"} or
    \code{"factors"}), the names by variable \code{label}. These scores
    are incompatible with conventional \code{plot} functions, but they can
    be used in \pkg{ggplot2}.}
    
  \item{labels}{Change plotting labels. The argument should be a list
    with elements \code{vectors} and \code{factors} which give the new
    plotting labels. If either of these elements is omitted, the
    default labels will be used. If there is only one type of elements
    (only \code{vectors} or only \code{factors}), the labels can be
    given as vector. The default labels can be displayed with
    \code{labels} command.}
  \item{arrow.mul}{Multiplier for vector lengths. The arrows are
    automatically scaled similarly as in \code{\link{plot.cca}} if this
    is not given in \code{plot} and \code{add = TRUE}. However, in
    \code{scores} it can be used to adjust arrow lengths when the
    \code{plot} function is not used.}
  \item{at}{The origin of fitted arrows in the plot.  If you plot arrows
    in other places then origin, you probably have to specify
    \code{arrrow.mul}.}
  \item{axis}{Plot axis showing the scaling of fitted arrows.}
  \item{p.max}{Maximum estimated \eqn{P} value for displayed
    variables.  You must calculate \eqn{P} values with setting
    \code{permutations} to use this option. }
  \item{col}{Colour in plotting.}
  \item{bg}{Background colour for labels. If \code{bg} is set, the
    labels are displayed with \code{\link{ordilabel}} instead of
    \code{text}. See Examples for using semitransparent background.}
  \item{add}{Results added to an existing ordination plot.}
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
  \item{display}{In fitting functions these are ordinary site scores or
    linear combination scores 
    (\code{"lc"}) in constrained ordination (\code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}}). In \code{scores}
    function they are either \code{"vectors"} or \code{"factors"}
    (with synonyms \code{"bp"} or \code{"cn"}, resp.).}
  \item{w}{Weights used in fitting (concerns mainly \code{\link{cca}}
    and \code{\link{decorana}} results which have nonconstant weights).}

  \item{...}{Parameters passed to \code{\link{scores}}.}
}
\details{
  Function \code{envfit} finds vectors or factor averages of
  environmental variables.  Function \code{plot.envfit} adds these in an
  ordination diagram.  If \code{X} is a \code{\link{data.frame}},
  \code{envfit}
  uses \code{factorfit} for \code{\link{factor}} variables and
  \code{vectorfit} for other variables.  If \code{X} is a matrix or a
  vector, \code{envfit} uses only \code{vectorfit}. Alternatively, the
  model can be defined a simplified model \code{\link{formula}}, where
  the left hand side must be an ordination result object or a matrix of
  ordination scores, and right hand
  side lists the environmental variables. The formula interface can be
  used for easier selection and/or transformation of environmental
  variables. Only the main effects will be analysed even if interaction
  terms were defined in the formula.

  The ordination results are extracted with \code{\link{scores}} and
  all extra arguments are passed to the \code{scores}. The fitted
  models only apply to the results defined when extracting the scores
  when using \code{envfit}. For instance, \code{scaling} in
  constrained ordination (see \code{\link{scores.rda}},
  \code{\link{scores.cca}}) must be set in the same way in
  \code{envfit} and in the \code{plot} or the ordination results (see
  Examples).

  The printed output of continuous variables (vectors) gives the
  direction cosines which are the coordinates of the heads of unit
  length vectors.  In \code{plot} these are scaled by their
  correlation (square root of the column \code{r2}) so that
  \dQuote{weak} predictors have shorter arrows than \dQuote{strong}
  predictors.  You can see the scaled relative lengths using command
  \code{scores}.  The \code{plot}ted (and scaled) arrows are further
  adjusted to the current graph using a constant multiplier: this will
  keep the relative \code{r2}-scaled lengths of the arrows but tries
  to fill the current plot.  You can see the multiplier using
  \code{ordiArrowMul(result_of_envfit)}, and set it with the
  argument \code{arrow.mul}.
  
  Functions \code{vectorfit} and \code{factorfit} can be called directly.
  Function \code{vectorfit} finds directions in the ordination space
  towards which the environmental vectors change most rapidly and to
  which they have maximal correlations with the ordination
  configuration.  Function \code{factorfit} finds averages of ordination
  scores for factor levels. Function \code{factorfit} treats ordered
  and unordered factors similarly.

  If \code{permutations} \eqn{> 0}, the significance of fitted vectors
  or factors is assessed using permutation of environmental variables.
  The goodness of fit statistic is squared correlation coefficient
  (\eqn{r^2}).
  For factors this is defined as \eqn{r^2 = 1 - ss_w/ss_t}, where
  \eqn{ss_w} and \eqn{ss_t} are within-group and total sums of
  squares. See \code{\link{permutations}} for additional details on
  permutation tests in Vegan.

  User can supply a vector of prior weights \code{w}. If the ordination
  object has weights, these will be used. In practise this means that
  the row totals are used as weights with \code{\link{cca}} or
  \code{\link{decorana}} results. If you do not like this, but want to
  give equal weights to all sites, you should set \code{w = NULL}.  The
  fitted vectors are similar to biplot arrows in constrained ordination
  only when fitted to LC scores (\code{display = "lc"}) and you set
  \code{scaling = "species"} (see \code{\link{scores.cca}}).  The
  weighted fitting gives similar results to biplot arrows and class
  centroids in \code{\link{cca}}.

  The lengths of arrows for fitted vectors are automatically adjusted
  for the physical size of the plot, and the arrow lengths cannot be
  compared across plots. For similar scaling of arrows, you must
  explicitly set the \code{arrow.mul} argument in the \code{plot}
  command; see \code{\link{ordiArrowMul}} and
  \code{\link{ordiArrowTextXY}}.
  
  The results can be accessed with \code{scores.envfit} function which
  returns either the fitted vectors scaled by correlation coefficient or
  the centroids of the fitted environmental variables, or a named list
  of both.
}
\value{
  Functions \code{vectorfit} and \code{factorfit} return lists of
  classes \code{vectorfit} and \code{factorfit} which have a
  \code{print} method.  The result object have the following items:

  \item{arrows}{Arrow endpoints from \code{vectorfit}. The arrows are
    scaled to unit length.}
  \item{centroids}{Class centroids from \code{factorfit}.}
  \item{r}{Goodness of fit statistic: Squared correlation coefficient}
  \item{permutations}{Number of permutations.}
  \item{control}{A list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}.}
  \item{pvals}{Empirical P-values for each variable.}

  Function \code{envfit} returns a list of class \code{envfit} with
  results of \code{vectorfit} and \code{envfit} as items.
  
  Function \code{plot.envfit} scales the vectors by correlation.
}

\author{Jari Oksanen.  The permutation test derives from the code
  suggested by Michael Scroggie. }

\note{
  Fitted vectors have become the method of choice in displaying
  environmental variables in ordination.  Indeed, they are the optimal
  way of presenting environmental variables in Constrained
  Correspondence Analysis \code{\link{cca}}, since there they are the
  linear constraints.
  In unconstrained ordination the relation between external variables
  and ordination configuration may be less linear, and therefore other
  methods than arrows may be more useful.  The simplest is to adjust the
  plotting symbol sizes (\code{cex}, \code{\link{symbols}}) by
  environmental variables.
  Fancier methods involve smoothing and regression methods that
  abound in \R, and \code{\link{ordisurf}} provides a wrapper for some.
  }

\seealso{
  A better alternative to vectors may be \code{\link{ordisurf}}.    
  }

\examples{
data(varespec, varechem)
library(MASS)
ord <- metaMDS(varespec)
(fit <- envfit(ord, varechem, perm = 999))
scores(fit, "vectors")
plot(ord)
plot(fit)
plot(fit, p.max = 0.05, col = "red")
## Adding fitted arrows to CCA. We use "lc" scores, and hope
## that arrows are scaled similarly in cca and envfit plots
ord <- cca(varespec ~ Al + P + K, varechem)
plot(ord, type="p")
fit <- envfit(ord, varechem, perm = 999, display = "lc")
plot(fit, p.max = 0.05, col = "red")
## 'scaling' must be set similarly in envfit and in ordination plot
plot(ord, type = "p", scaling = "sites")
fit <- envfit(ord, varechem, perm = 0, display = "lc", scaling = "sites")
plot(fit, col = "red")

## Class variables, formula interface, and displaying the
## inter-class variability with ordispider, and semitransparent
## white background for labels (semitransparent colours are not
## supported by all graphics devices)
data(dune)
data(dune.env)
ord <- cca(dune)
fit <- envfit(ord ~ Moisture + A1, dune.env, perm = 0)
plot(ord, type = "n")
with(dune.env, ordispider(ord, Moisture, col="skyblue"))
with(dune.env, points(ord, display = "sites", col = as.numeric(Moisture),
                      pch=16))
plot(fit, cex=1.2, axis=TRUE, bg = rgb(1, 1, 1, 0.5))
## Use shorter labels for factor centroids
labels(fit)
plot(ord)
plot(fit, labels=list(factors = paste("M", c(1,2,4,5), sep = "")),
   bg = rgb(1,1,0,0.5))
}
\keyword{multivariate }
\keyword{aplot}
\keyword{htest}




