\name{capscale}
\alias{capscale}
\alias{print.capscale}

\title{[Partial] Constrained Analysis of Principal Coordinates or
  distance-based RDA }
\description{
  Constrained Analysis of Principal Coordinates (CAP) is an ordination method
  similar to Redundancy Analysis (\code{\link{rda}}), but it allows
  non-Euclidean dissimilarity indices, such as Manhattan or
  Bray--Curtis distance. Despite this non-Euclidean feature, the analysis
  is strictly linear and metric. If called with Euclidean distance,
  the results are identical to \code{\link{rda}}, but \code{capscale}
  will be much more inefficient. Function \code{capscale} is a
  constrained version of metric scaling, a.k.a. principal coordinates
  analysis, which is based on the Euclidean distance but can be used,
  and is more useful, with other dissimilarity measures. The function
  can also perform unconstrained principal coordinates analysis,
  optionally using extended dissimilarities.
}
\usage{
capscale(formula, data, distance = "euclidean", sqrt.dist = FALSE,
    comm = NULL, add = FALSE,  dfun = vegdist, metaMDSdist = FALSE,
    na.action = na.fail, subset = NULL, ...)
}

\arguments{
  \item{formula}{Model formula. The function can be called only with the
    formula interface. Most usual features of \code{\link{formula}} hold,
    especially as defined in \code{\link{cca}} and \code{\link{rda}}. The
    LHS must be either a community data matrix or a dissimilarity matrix,
    e.g., from
    \code{\link{vegdist}} or \code{\link{dist}}.
    If the LHS is a data matrix, function \code{\link{vegdist}}
    will be used to find the dissimilarities. The RHS defines the
    constraints. The constraints can be continuous variables or factors,
    they can be transformed within the formula, and they can have
    interactions as in a typical \code{\link{formula}}. The RHS can have a
    special term \code{Condition} that defines variables to be
    ``partialled out'' before constraints, just like in \code{\link{rda}}
    or \code{\link{cca}}. This allows the use of partial CAP.}
  \item{data}{ Data frame containing the variables on the right hand side of
    the model formula. }
  \item{distance}{The name of the dissimilarity (or distance) index if
    the LHS of the \code{formula} is a data frame instead of
    dissimilarity matrix.}
  \item{sqrt.dist}{Take square roots of dissimilarities. See section
  \code{Notes} below.}
  \item{comm}{ Community data frame which will be used for finding
    species scores when the LHS of the \code{formula} was a
    dissimilarity matrix. This is not used if the LHS is a data
    frame. If this is not supplied, the ``species scores'' are the axes
    of initial metric scaling (\code{\link{cmdscale}}) and may be
    confusing.}
  \item{add}{Logical indicating if an additive constant should be
     computed, and added to the non-diagonal dissimilarities such
     that all eigenvalues are non-negative in the underlying
     Principal Co-ordinates Analysis (see \code{\link{cmdscale}} 
     for details). This implements \dQuote{correction method 2} of
     Legendre & Legendre (1998, p. 434). The negative eigenvalues are
     caused by using semi-metric or non-metric dissimilarities with
     basically metric \code{\link{cmdscale}}. They are harmless and
     ignored in \code{capscale}, but you also can avoid warnings with
     this option. }
   \item{dfun}{Distance or dissimilarity function used. Any function
     returning standard \code{"dist"} and taking the index name as the
     first argument can be used. }
   \item{metaMDSdist}{Use \code{\link{metaMDSdist}} similarly as in
     \code{\link{metaMDS}}. This means automatic data transformation and
     using extended flexible shortest path dissimilarities (function
     \code{\link{stepacross}}) when there are many dissimilarities based on
     no shared species.}
   \item{na.action}{Handling of missing values in constraints or
     conditions. The default (\code{\link{na.fail}}) is to stop
     with missing values. Choices \code{\link{na.omit}} and
     \code{\link{na.exclude}} delete rows with missing values, but
     differ in representation of results. With \code{na.omit} only
     non-missing site scores are shown, but \code{na.exclude} gives
     \code{NA} for scores of missing observations. Unlike in
     \code{\link{rda}}, no WA scores are available for missing
     constraints or conditions.}
   \item{subset}{Subset of data rows. This can be a logical vector
     which is \code{TRUE} for kept observations, or a logical
     expression which can contain variables in the working
     environment, \code{data} or species names of the community data
     (if given in the formula or as \code{comm} argument).}
   \item{\dots}{Other parameters passed to \code{\link{rda}} or to
     \code{\link{metaMDSdist}}.  }
}
\details{
  Canonical Analysis of Principal Coordinates (CAP) is simply a
  Redundancy Analysis of results of Metric (Classical) Multidimensional
  Scaling (Anderson & Willis 2003). Function capscale uses two steps:
  (1) it ordinates the dissimilarity matrix using
  \code{\link{cmdscale}} and (2) analyses these results using
  \code{\link{rda}}. If the user supplied a community data frame instead
  of dissimilarities, the function will find the needed dissimilarity
  matrix using \code{\link{vegdist}} with specified
  \code{distance}. However, the method will accept dissimilarity
  matrices from \code{\link{vegdist}}, \code{\link{dist}}, or any
  other method producing similar matrices. The constraining variables can be
  continuous or factors or both, they can have interaction terms,
  or they can be transformed in the call. Moreover, there can be a
  special term
  \code{Condition} just like in \code{\link{rda}} and \code{\link{cca}}
  so that ``partial'' CAP can be performed.

  The current implementation  differs from the method suggested by
  Anderson & Willis (2003) in three major points which actually make it
  similar to distance-based redundancy analysis (Legendre & Anderson
  1999):
  \enumerate{
    \item Anderson & Willis used the orthonormal solution of
    \code{\link{cmdscale}}, whereas \code{capscale} uses axes
    weighted by corresponding eigenvalues, so that the ordination
    distances are the best approximations of original
    dissimilarities. In the original method, later ``noise'' axes are
    just as important as first major axes.
    \item Anderson & Willis take only a subset of axes, whereas 
    \code{capscale} uses all axes with positive eigenvalues. The use of
    subset is necessary with orthonormal axes to chop off some
    ``noise'', but the use of all axes guarantees that the results are
    the best approximation of original dissimilarities.
    \item Function \code{capscale} adds species scores as weighted sums
    of (residual) community matrix (if the matrix is available), whereas
    Anderson & Willis have no fixed method for adding species scores.
  }
  With these definitions, function \code{capscale} with Euclidean
  distances will be identical to \code{\link{rda}} in eigenvalues and
  in site, species and biplot scores (except for possible sign
  reversal). 
  However, it makes no sense to use \code{capscale} with
  Euclidean distances, since direct use of \code{\link{rda}} is much more
  efficient. Even with non-Euclidean dissimilarities, the
  rest of the analysis will be metric and linear.

  The function can be also used to perform ordinary metric scaling
  a.k.a. principal coordinates analysis by using a formula with only a
  constant on the left hand side, or \code{comm ~ 1}. With
  \code{metaMDSdist = TRUE}, the function can do automatic data
  standardization and use extended dissimilarities using function
  \code{\link{stepacross}} similarly as in non-metric multidimensional
  scaling with \code{\link{metaMDS}}. 
  
}
\value{
  The function returns an object of class \code{capscale} which is
  identical to the result of \code{\link{rda}}. At the moment,
  \code{capscale} does not have specific methods, but it uses
  \code{\link{cca}} and \code{\link{rda}} methods
  \code{\link{plot.cca}},
  \code{\link{scores.rda}}  etc. Moreover, you
  can use \code{\link{anova.cca}} for permutation tests of
  ``significance'' of the results.
}
\references{
  Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal
  coordinates: a useful method of constrained ordination for
  ecology. \emph{Ecology} 84, 511--525.

  Gower, J.C. (1985). Properties of Euclidean and non-Euclidean
  distance matrices. \emph{Linear Algebra and its Applications} 67, 81--97.

  Legendre, P. & Anderson, M. J. (1999). Distance-based redundancy
  analysis: testing multispecies responses in multifactorial ecological
  experiments. \emph{Ecological Monographs} 69, 1--24.

  Legendre, P. & Legendre, L. (1998).  \emph{Numerical Ecology}. 2nd English
  Edition. Elsevier
}
\author{ Jari Oksanen }

\note{ The function produces negative eigenvalues with non-Euclidean
  dissimilarity indices. The non-Euclidean component of inertia is
  given under the title \code{Imaginary} in the printed output. The
  \code{Total} inertia is the sum of all eigenvalues, but the sum of
  all non-negative eigenvalues is given as \code{Real Total} (which is
  higher than the \code{Total}). The ordination is based only on the
  real dimensions with positive eigenvalues, and therefore the
  proportions of inertia components only apply to the \code{Real
  Total} and ignore the \code{Imaginary} component. Permutation tests
  with \code{\link{anova.cca}} use only the real solution of positive
  eigenvalues. Function \code{\link{adonis}} gives similar
  significance tests, but it also handles the imaginary dimensions
  (negative eigenvalues) and therefore its results may differ from
  permutation test results of \code{capscale}.

  If the negative eigenvalues are disturbing, you can
  use argument \code{add = TRUE} passed to \code{\link{cmdscale}}, or,
  preferably, a distance measure that does not cause these warnings.
  Alternatively, after square root transformation of distances
  (argument \code{sqrt.dist = TRUE}) many indices do not produce
  negative eigenvalues.

  The inertia is named after the dissimilarity index as defined in the
  dissimilarity data, or as \code{unknown distance} if such an
  information is missing.  Function \code{\link{rda}} usually divides
  the ordination scores by number of sites minus one. In this way, the
  inertia is variance instead of sum of squares, and the eigenvalues sum
  up to variance. Many dissimilarity measures are in the range 0 to 1,
  so they have already made a similar division. If the largest original
  dissimilarity is less than or equal to 4 (allowing for
  \code{\link{stepacross}}), this division is undone in \code{capscale}
  and original dissimilarities are used. Keyword \code{mean} is added to
  the inertia in cases where division was made, e.g. in Euclidean and
  Manhattan distances.  Inertia is based on squared index, and keyword
  \code{squared} is added to the name of distance, unless data were
  square root transformed (argument \code{sqrt.dist = TRUE}). If an
  additive constant was used, keyword \code{euclidified} is added to the
  the name of inertia (argument \code{add = TRUE}).
}


\seealso{\code{\link{rda}}, \code{\link{cca}}, \code{\link{plot.cca}},
  \code{\link{anova.cca}}, \code{\link{vegdist}},
  \code{\link{dist}}, \code{\link{cmdscale}}.}
\examples{
data(varespec)
data(varechem)
## Basic Analysis
vare.cap <- capscale(varespec ~ N + P + K + Condition(Al), varechem,
                     dist="bray")
vare.cap
plot(vare.cap)
anova(vare.cap)
## Avoid negative eigenvalues with additive constant
capscale(varespec ~ N + P + K + Condition(Al), varechem,
                     dist="bray", add =TRUE)
## Avoid negative eigenvalues by taking square roots of dissimilarities
capscale(varespec ~ N + P + K + Condition(Al), varechem,
                     dist = "bray", sqrt.dist= TRUE)
## Principal coordinates analysis with extended dissimilarities
capscale(varespec ~ 1, dist="bray", metaMDS = TRUE)
}
\keyword{ multivariate }

