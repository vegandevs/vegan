\name{anosim}
\alias{anosim}
\alias{print.anosim}
\alias{summary.anosim}
\alias{plot.anosim}

\title{ Analysis of Similarities }
\description{
  Analysis of similarities (ANOSIM) provides a way to test statistically
  whether there is a significant difference between two or more groups
  of sampling units.
}
\usage{
anosim(dat, grouping, permutations = 999, distance = "bray", strata)
}

\arguments{
  \item{dat}{Data matrix or data frame in which rows are samples and
    columns are response variable(s), or a dissimilarity object or a
    symmetric square matrix of dissimilarities.}
  \item{grouping}{Factor for grouping observations.}
  \item{permutations}{Number of permutation to assess the significance
    of the ANOSIM statistic. }
  \item{distance}{Choice of distance metric that measures the
    dissimilarity between two observations . See \code{\link{vegdist}} for
    options.  This will be used if \code{dat} was not a dissimilarity
    structure or a symmetric square matrix.}  
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
}
\details{
  Analysis of similarities (ANOSIM) provides a way to test statistically
  whether there is a significant difference between two or more groups
  of sampling units.  Function \code{anosim} operates directly on a
  dissimilarity matrix.  A suitable dissimilarity matrix is produced by
  functions \code{\link{dist}} or \code{\link{vegdist}}.  The
  method is philosophically allied with NMDS ordination
  (\code{\link[MASS]{isoMDS}}), in that it uses only the rank order of
  dissimilarity values.

  If two groups of sampling units are really different in their species
  composition, then compositional dissimilarities between the groups
  ought to be greater than those within the groups.  The \code{anosim}
  statistic \eqn{R} is based on the difference of mean ranks between
  groups (\eqn{r_B}) and within groups (\eqn{r_W}):

  \deqn{R = (r_B - r_W)/(N (N-1) / 4)}

  The divisor is chosen so that \eqn{R} will be in the interval
  \eqn{-1 \dots +1}, value \eqn{0} indicating completely random
  grouping.

  The statistical significance of observed \eqn{R} is assessed by
  permuting the grouping vector to obtain the empirical
  distribution of \eqn{R} under null-model.

  The function has \code{summary} and \code{plot} methods.  These both
  show valuable information to assess the validity of the method:  The
  function assumes that all ranked dissimilarities within groups 
  have about equal median and range.  The \code{plot} method uses
  \code{\link{boxplot}} with options \code{notch=TRUE} and
  \code{varwidth=TRUE}. 
}
\value{
  The function returns a list of class \code{"anosim"} with following
  items: 
  \item{call }{Function call.}
  \item{statistic}{The value of ANOSIM statistic \eqn{R}}
  \item{signif}{Significance from permutation.}
  \item{perm}{Permutation values of \eqn{R}}
  \item{class.vec}{Factor with value \code{Between} for dissimilarities
    between classes and class name for corresponding dissimilarity
    within class.}
  \item{dis.rank}{Rank of dissimilarity entry.}
  \item{dissimilarity}{The name of the dissimilarity index: the
    \code{"method"} entry of the \code{dist} object.}
}
\references{
  Clarke, K. R. (1993). Non-parametric multivariate analysis of changes
  in community structure. \emph{Australian Journal of Ecology} 18, 117-143.
}
\author{Jari Oksanen, with a help from Peter R. Minchin.}
\note{
  I don't quite trust this method.  Somebody should study its
  performance carefully.  The function returns a lot of information 
  to ease further scrutiny. Most \code{anosim} models could be analysed
  with \code{\link{adonis}} which seems to be a more robust alternative.

}

\seealso{\code{\link{mrpp}} for a similar function using original
dissimilarities instead of their ranks. 
\code{\link{dist}} and \code{\link{vegdist}} for obtaining
  dissimilarities, and \code{\link{rank}} for ranking real values.  For
  comparing dissimilarities against continuous variables, see
  \code{\link{mantel}}. Function \code{\link{adonis}} is a more robust
  alternative that should preferred. }

\examples{
data(dune)
data(dune.env)
dune.dist <- vegdist(dune)
attach(dune.env)
dune.ano <- anosim(dune.dist, Management)
summary(dune.ano)
plot(dune.ano)
}
\keyword{multivariate }
\keyword{ nonparametric }
\keyword{ htest }
