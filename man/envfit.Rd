\name{envfit}
\alias{envfit}
\alias{envfit.default}
\alias{envfit.formula}
\alias{vectorfit}
\alias{factorfit}
\alias{plot.envfit}
\alias{print.envfit}
\alias{print.factorfit}
\alias{print.vectorfit}
\alias{scores.envfit}

\title{Fits an Environmental Vector or Factor onto an Ordination }
\description{
  The function fits environmental vectors or factors onto an
  ordination. The projections of points onto vectors have
  maximum correlation with corresponding environmental variables, and
  the factors show the averages of factor levels.
}
\usage{
\method{envfit}{default}(ord, env, permutations = 999, strata, choices=c(1,2), 
   display = "sites", w  = weights(ord), na.rm = FALSE, ...)
\method{envfit}{formula}(formula, data, ...)
\method{plot}{envfit}(x, choices = c(1,2), arrow.mul, at = c(0,0), axis = FALSE, 
    p.max = NULL, col = "blue", add = TRUE, ...)
\method{scores}{envfit}(x, display, choices, ...)
vectorfit(X, P, permutations = 0, strata, w, ...)
factorfit(X, P, permutations = 0, strata, w, ...)
}

\arguments{
  \item{ord}{An ordination object or other structure from which the
    ordination \code{\link{scores}} can be extracted (including a data
    frame or matrix of scores)}.
  \item{env}{Data frame, matrix or vector of environmental
    variables. The variables can be of mixed type (factors, continuous
    variables) in data frames.}
  \item{X}{Matrix or data frame of ordination scores.}
  \item{P}{Data frame, matrix or vector of environmental
    variable(s). These must be continuous for \code{vectorfit} and
    factors or characters for \code{factorfit}. }
  \item{permutations}{ Number of permutations for assessing significance
    of vectors or factors. Set to \code{0} to skip permutations.}
  \item{formula, data}{Model  \code{\link{formula}} and data.  }
  \item{na.rm}{Remove points with missing values in ordination scores
    or environmental variables. The operation is casewise: the whole
    row of data is removed if there is a missing value and 
    \code{na.rm = TRUE}.}
  \item{x}{A result object from \code{envfit}.}
  \item{choices}{Axes to plotted.}
  \item{arrow.mul}{Multiplier for vector lengths. The arrows are
    automatically scaled similarly as in \code{\link{plot.cca}} if this
    is not given and \code{add = TRUE}.}
  \item{at}{The origin of fitted arrows in the plot.  If you plot arrows
    in other places then origin, you probably have to specify
    \code{arrrow.mul}.}
  \item{axis}{Plot axis showing the scaling of fitted arrows.}
  \item{p.max}{Maximum estimated \eqn{P} value for displayed
    variables.  You must calculate \eqn{P} values with setting
    \code{permutations} to use this option. }
  \item{col}{Colour in plotting.}
  \item{add}{Results added to an existing ordination plot.}
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
  \item{display}{In fitting functions these are ordinary site scores or
    linear combination scores 
    (\code{"lc"}) in constrained ordination (\code{\link{cca}},
    \code{\link{rda}}, \code{\link{capscale}}). In \code{scores}
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
  
  Functions \code{vectorfit} and \code{factorfit} can be called directly.
  Function \code{vectorfit} finds directions in the ordination space
  towards which the environmental vectors change most rapidly and to
  which they have maximal correlations with the ordination
  configuration.  Function \code{factorfit} finds averages of ordination
  scores for factor levels. Function \code{factorfit} treats ordered
  and unordered factors similarly.

  If \code{permutations} \eqn{> 0}, the `significance' of fitted vectors
  or factors is assessed using permutation of environmental variables.
  The goodness of fit statistic is squared correlation coefficient
  (\eqn{r^2}).
  For factors this is defined as \eqn{r^2 = 1 - ss_w/ss_t}, where
  \eqn{ss_w} and \eqn{ss_t} are within-group and total sums of
  squares. See \code{\link{permutations}} for additional details on
  permutation tests in Vegan.

  User can supply a vector of prior  weights \code{w}. If the ordination
  object has weights, these will be used. In practise this means that
  the row totals are used as weights with
  \code{\link{cca}} or \code{\link{decorana}} results. If you do not
  like this, but want to give 
  equal weights to all sites, you should set \code{w = NULL}.
  The weighted fitting gives similar results to biplot
  arrows and class centroids in \code{\link{cca}}.
  For complete
  similarity between fitted vectors and biplot arrows, you should set
  \code{display = "lc"} (and possibly \code{scaling = 2}).

  The lengths of arrows for fitted vectors are automatically adjusted
  for the physical size of the plot, and the arrow lengths cannot be
  compared across plots. For similar scaling of arrows, you must
  explicitly set the \code{arrow.mul} argument in the \code{plot}
  command.

  The results can be accessed with \code{scores.envfit} function which
  returns either the fitted vectors scaled by correlation coefficient or
  the centroids of the fitted environmental variables.
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
data(varespec)
data(varechem)
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
## Class variables, formula interface, and displaying the
## inter-class variability with `ordispider'
data(dune)
data(dune.env)
attach(dune.env)
ord <- cca(dune)
fit <- envfit(ord ~ Moisture + A1, dune.env)
plot(ord, type = "n")
ordispider(ord, Moisture, col="skyblue")
points(ord, display = "sites", col = as.numeric(Moisture), pch=16)
plot(fit, cex=1.2, axis=TRUE)
}
\keyword{multivariate }
\keyword{aplot}
\keyword{htest}




