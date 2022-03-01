\name{stepacross}
\alias{stepacross}

\title{Stepacross as Flexible Shortest Paths or Extended Dissimilarities } 
\description{
  Function \code{stepacross} tries to replace dissimilarities with
  shortest paths stepping across intermediate 
  sites while regarding dissimilarities above a threshold as missing
  data (\code{NA}). With \code{path = "shortest"} this is the flexible shortest
  path (Williamson 1978, Bradfield & Kenkel 1987),
  and with \code{path = "extended"} an
  approximation known as extended dissimilarities (De'ath 1999).
  The use of \code{stepacross} should improve the ordination with high
  beta diversity, when there are many sites with no species in common.
}
\usage{
stepacross(dis, path = "shortest", toolong = 1, trace = TRUE, ...)
}
\arguments{
  \item{dis}{Dissimilarity data inheriting from class \code{dist} or
    a an object, such as a matrix, that can be converted to a
    dissimilarity matrix. Functions \code{\link{vegdist}} and
    \code{\link{dist}} are some functions producing suitable
    dissimilarity data. }
  \item{path}{The method of stepping across (partial match)
    Alternative \code{"shortest"} finds the shortest paths, and
    \code{"extended"}  their approximation known as extended
    dissimilarities.} 
  \item{toolong}{Shortest dissimilarity regarded as \code{NA}.
    The function uses a fuzz factor, so
    that dissimilarities close to the limit will be made \code{NA}, too. }
  \item{trace}{ Trace the calculations.}
  \item{\dots}{Other parameters (ignored).}
}
\details{
  Williamson (1978) suggested using flexible shortest paths to estimate
  dissimilarities between sites which have nothing in common, or no shared
  species. With \code{path = "shortest"} function \code{stepacross}
  replaces dissimilarities that are
  \code{toolong} or longer with \code{NA}, and tries to find shortest
  paths between all sites using remaining dissimilarities. Several
  dissimilarity indices are semi-metric which means that they do not
  obey the triangle inequality \eqn{d_{ij} \leq d_{ik} + d_{kj}}{d[ij] <=
    d[ik] + d[kj]}, and shortest path algorithm can replace these
  dissimilarities as well, even when they are shorter than
  \code{toolong}. 

  De'ath (1999) suggested a simplified method known as extended
  dissimilarities, which are calculated with \code{path = "extended"}. 
  In this method, dissimilarities that are
  \code{toolong} or longer are first made \code{NA}, and then the function
  tries to replace these \code{NA} dissimilarities with a path through
  single stepping stone points. If not all \code{NA} could be 
  replaced with one pass, the function will make new passes with updated
  dissimilarities as long as
  all \code{NA} are replaced with extended dissimilarities. This mean
  that in the second and further passes, the remaining \code{NA}
  dissimilarities are allowed to have more than one stepping stone site,
  but previously replaced dissimilarities are not updated. Further, the
  function does not consider dissimilarities shorter than \code{toolong},
  although some of these could be replaced with a shorter path in
  semi-metric indices, and used as a part of other paths. In optimal
  cases, the extended dissimilarities are equal to shortest paths, but
  they may be longer.  

  As an alternative to defining too long dissimilarities with parameter
  \code{toolong}, the input dissimilarities can contain \code{NA}s. If
  \code{toolong} is zero or negative, the function does not make any
  dissimilarities into \code{NA}. If there are no \code{NA}s in the
  input  and \code{toolong = 0}, \code{path = "shortest"}
  will find shorter paths for semi-metric indices, and \code{path = "extended"} 
  will do nothing. Function \code{\link{no.shared}} can be
  used to set dissimilarities to \code{NA}.
  
  If the data are disconnected or there is no path between all points,
  the result will
  contain \code{NA}s and a warning is issued. Several methods cannot
  handle \code{NA} dissimilarities, and this warning should be taken
  seriously. Function \code{\link{distconnected}} can be used to find
  connected groups and remove rare outlier observations or groups of
  observations.

  Alternative \code{path = "shortest"} uses Dijkstra's method for
  finding flexible shortest paths, implemented as priority-first search
  for dense graphs (Sedgewick 1990). Alternative \code{path = "extended"} 
  follows De'ath (1999), but implementation is simpler
  than in his code.
  
}
\value{
  Function returns an object of class \code{dist} with extended
  dissimilarities (see functions \code{\link{vegdist}} and
  \code{\link{dist}}). 
  The value of \code{path} is appended to the \code{method} attribute.
}
\references{
  Bradfield, G.E. & Kenkel, N.C. (1987). Nonlinear ordination using
  flexible shortest path adjustment of ecological
  distances. \emph{Ecology} 68, 750--753.
  
  De'ath, G. (1999). Extended dissimilarity: a method of robust
  estimation of ecological distances from high beta diversity
  data. \emph{Plant Ecol.} 144, 191--199.

  Sedgewick, R. (1990). \emph{Algorithms in C}. Addison Wesley. 

  Williamson, M.H. (1978). The ordination of incidence
  data. \emph{J. Ecol.} 66, 911-920.
}
\author{ Jari Oksanen}
\note{
  The function changes the original dissimilarities, and not all
  like this. It may be best to  use  the
  function only when you really \emph{must}:  extremely high
  beta diversity where a large proportion of dissimilarities are at their
  upper limit (no species in common). 

  Semi-metric indices vary in their degree of violating the triangle
  inequality. Morisita and Horn--Morisita indices of
  \code{\link{vegdist}} may be very strongly semi-metric, and shortest
  paths can change these indices very much. Mountford index violates
  basic rules of dissimilarities: non-identical sites have zero
  dissimilarity if species composition of the poorer site is a subset of
  the richer. With Mountford index, you can find three sites \eqn{i, j,
    k} so that \eqn{d_{ik} = 0}{d[ik] = 0} and \eqn{d_{jk} = 0}{d[jk] =
    0}, but \eqn{d_{ij} > 0}{d[ij] > 0}. The results of \code{stepacross}
  on Mountford index can be very weird. If \code{stepacross} is needed,
  it is best to try to use it with more metric indices only.
}

\seealso{
  Function \code{\link{distconnected}} can find connected groups in
  disconnected data, and function \code{\link{no.shared}} can be used to
  set dissimilarities as \code{NA}.  See \code{\link{swan}} for an
  alternative approach. Function \code{stepacross} is an essential
  component in \code{\link{isomap}} and \code{\link{cophenetic.spantree}}.
 }
\examples{
# There are no data sets with high beta diversity in vegan, but this
# should give an idea.
data(dune)
dis <- vegdist(dune)
edis <- stepacross(dis)
plot(edis, dis, xlab = "Shortest path", ylab = "Original")
## Manhattan distance have no fixed upper limit.
dis <- vegdist(dune, "manhattan")
is.na(dis) <- no.shared(dune)
dis <- stepacross(dis, toolong=0)
}
\keyword{multivariate }

