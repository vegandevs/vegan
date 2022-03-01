\name{capscale}
\alias{capscale}
\alias{dbrda}

\title{[Partial] Distance-based Redundancy Analysis }
\description{
  Distance-based redundancy analysis (dbRDA) is an ordination method
  similar to Redundancy Analysis (\code{\link{rda}}), but it allows
  non-Euclidean dissimilarity indices, such as Manhattan or
  Bray--Curtis distance. Despite this non-Euclidean feature, the analysis
  is strictly linear and metric. If called with Euclidean distance,
  the results are identical to \code{\link{rda}}, but dbRDA
  will be less efficient. Functions \code{capscale} and \code{dbrda} are
  constrained versions of metric scaling, a.k.a. principal coordinates
  analysis, which are based on the Euclidean distance but can be used,
  and are more useful, with other dissimilarity measures. The functions
  can also perform unconstrained principal coordinates analysis,
  optionally using extended dissimilarities.
}
\usage{
capscale(formula, data, distance = "euclidean", sqrt.dist = FALSE,
    comm = NULL, add = FALSE,  dfun = vegdist, metaMDSdist = FALSE,
    na.action = na.fail, subset = NULL, ...)
dbrda(formula, data, distance = "euclidean", sqrt.dist = FALSE,
    add = FALSE, dfun = vegdist, metaMDSdist = FALSE,
    na.action = na.fail, subset = NULL, ...)
}

\arguments{
  \item{formula}{Model formula. The function can be called only with the
    formula interface. Most usual features of \code{\link{formula}} hold,
    especially as defined in \code{\link{cca}} and \code{\link{rda}}. The
    LHS must be either a community data matrix or a dissimilarity matrix,
    e.g., from
    \code{\link{vegdist}} or \code{\link{dist}}.
    If the LHS is a data matrix, function \code{\link{vegdist}} or
    function given in \code{dfun}
    will be used to find the dissimilarities. The RHS defines the
    constraints. The constraints can be continuous variables or factors,
    they can be transformed within the formula, and they can have
    interactions as in a typical \code{\link{formula}}. The RHS can have a
    special term \code{Condition} that defines variables to be
    \dQuote{partialled out} before constraints, just like in \code{\link{rda}}
    or \code{\link{cca}}. This allows the use of partial dbRDA.}
  \item{data}{ Data frame containing the variables on the right hand side of
    the model formula. }
  \item{distance}{The name of the dissimilarity (or distance) index if
    the LHS of the \code{formula} is a data frame instead of
    dissimilarity matrix.}
  \item{sqrt.dist}{Take square roots of dissimilarities. See section
    \code{Details} below.}

  \item{comm}{ Community data frame which will be used for finding
    species scores when the LHS of the \code{formula} was a
    dissimilarity matrix. This is not used if the LHS is a data
    frame. If this is not supplied, the \dQuote{species scores} are
    unavailable when dissimilarities were supplied. N.B., this is
    only available in \code{capscale}: \code{dbrda} does not return
    species scores. Function \code{\link{sppscores}} can be used to add
    species scores if they are missing.}

  \item{add}{Add a constant to the non-diagonal dissimilarities such
     that all eigenvalues are non-negative in the underlying Principal
     Co-ordinates Analysis (see \code{\link{wcmdscale}} for
     details). \code{"lingoes"} (or \code{TRUE}) uses the
     recommended method of Legendre & Anderson (1999: \dQuote{method
     1}) and \code{"cailliez"} uses their \dQuote{method 2}. The
     latter is the only one in \code{\link{cmdscale}}.}
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
   \item{\dots}{Other parameters passed to underlying functions (e.g.,
     \code{\link{metaMDSdist}}).  }
}
\details{

  Functions \code{capscale} and \code{dbrda} provide two alternative
  implementations of dbRDA. Function \code{capscale} is based on
  Legendre & Anderson (1999): the dissimilarity data are first
  ordinated using metric scaling, and the ordination results are
  analysed as \code{\link{rda}}. Function \code{dbrda} is based on
  McArdle & Anderson (2001) and directly decomposes
  dissimilarities. It does not use \code{\link{rda}} but a parallel
  implementation adapted for analysing dissimilarities and returns a
  subset of \code{\link{rda}} items. With Euclidean distances both
  results are identical to \code{\link{rda}}.  Other dissimilarities
  may give negative eigenvalues associated with imaginary
  axes. Negative eigenvalues are handled differently: \code{capscale}
  ignores imaginary axes and analyses only real axes with positive
  eigenvalues, and \code{dbrda} directly analyses dissimilarities and
  can give negative eigenvalues in any component.

  If the user supplied a community data frame instead of
  dissimilarities, the functions will find dissimilarities using
  \code{\link{vegdist}} or distance function given in \code{dfun} with
  specified \code{distance}. The functions will accept distance
  objects from \code{\link{vegdist}}, \code{\link{dist}}, or any other
  method producing compatible objects. The constraining variables can be
  continuous or factors or both, they can have interaction terms, or
  they can be transformed in the call. Moreover, there can be a
  special term \code{Condition} just like in \code{\link{rda}} and
  \code{\link{cca}} so that \dQuote{partial} analysis can be performed.

  Function \code{dbrda} does not return species scores, and they can
  also be missing in \code{capscale}, but they can be added after the
  analysis using function \code{\link{sppscores}}.

  Non-Euclidean dissimilarities can produce negative eigenvalues
  (Legendre & Anderson 1999, McArdle & Anderson 2001).  If there are
  negative eigenvalues, the printed output of \code{capscale} will add
  a column with sums of positive eigenvalues and an item of sum of
  negative eigenvalues, and \code{dbrda} will add a column giving the
  number of real dimensions with positive eigenvalues.  If negative
  eigenvalues are disturbing, functions let you to distort the
  dissimilarities so that only non-negative eigenvalues will be
  produced with argument \code{add = TRUE}. Alternatively, with
  \code{sqrt.dist = TRUE}, square roots of dissimilarities will be
  used which may help in avoiding negative eigenvalues (Legendre &
  Anderson 1999).

  The functions can be also used to perform ordinary metric scaling
  a.k.a. principal coordinates analysis by using a formula with only a
  constant on the left hand side, or \code{comm ~ 1}. With
  \code{metaMDSdist = TRUE}, the function can do automatic data
  standardization and use extended dissimilarities using function
  \code{\link{stepacross}} similarly as in non-metric multidimensional
  scaling with \code{\link{metaMDS}}.

}
\value{

  The functions return an object of class \code{capscale} or
  \code{dbrda} which inherits from \code{\link{rda}}. See
  \code{\link{cca.object}} for description of the result object.

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

  Legendre, P. & Legendre, L. (2012).  \emph{Numerical Ecology}. 3rd English
  Edition. Elsevier.

  McArdle, B.H. & Anderson, M.J. (2001). Fitting multivariate models
  to community data: a comment on distance-based redundancy
  analysis. \emph{Ecology} 82, 290--297.
}

\author{ Jari Oksanen }

\note{ The function \code{capscale} was originally developed as a
  variant of constrained analysis of proximities (Anderson & Willis
  2003), but these developments made it similar to dbRDA. However, it
  discards the imaginary dimensions with negative eigenvalues and
  ordination and significance tests area only based on real dimensions
  and positive eigenvalues.

  The inertia is named after the dissimilarity index as defined in the
  dissimilarity data, or as \code{unknown distance} if such
  information is missing. If the largest original dissimilarity was
  larger than 4, \code{capscale} handles input similarly as \code{rda}
  and bases its analysis on variance instead of sum of
  squares. Keyword \code{mean} is added to the inertia in these cases,
  e.g. with Euclidean and Manhattan distances.  Inertia is based on
  squared index, and keyword \code{squared} is added to the name of
  distance, unless data were square root transformed (argument
  \code{sqrt.dist=TRUE}). If an additive constant was used with
  argument \code{add}, \code{Lingoes} or \code{Cailliez adjusted} is
  added to the the name of inertia, and the value of the constant is
  printed.}


\seealso{\code{\link{rda}}, \code{\link{cca}}, \code{\link{plot.cca}},
  \code{\link{anova.cca}}, \code{\link{vegdist}},
  \code{\link{dist}}, \code{\link{cmdscale}}, \code{\link{wcmdscale}}
  for underlying and related functions. Function \code{\link{sppscores}}
  can add species scores or replace existing species scores.

  The function returns similar result object as \code{\link{rda}} (see
  \code{\link{cca.object}}). This section for \code{\link{rda}} gives a
  more complete list of functions that can be used to access and
  analyse dbRDA results.

}
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
## dbrda
dbrda(varespec ~ N + P + K + Condition(Al), varechem,
                     dist="bray")
## avoid negative eigenvalues also with Jaccard distances
dbrda(varespec ~ N + P + K + Condition(Al), varechem,
                     dist="jaccard")
}
\keyword{ multivariate }

