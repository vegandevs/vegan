\name{distconnected}
\alias{distconnected}
\alias{no.shared}

\title{Connectedness of Dissimilarities}
\description{
  Function \code{distconnected} finds groups that are connected
  disregarding dissimilarities that are at or above a threshold or
  \code{NA}. The function can be used to find groups that can be
  ordinated together or transformed by
  \code{\link{stepacross}}. Function \code{no.shared} returns a logical
  dissimilarity object, where \code{TRUE} means that sites have no
  species in common. This is a minimal structure for
  \code{distconnected} or can be used to set missing values to
  dissimilarities.
}
\usage{
distconnected(dis, toolong = 1, trace = TRUE)

no.shared(x)
}

\arguments{
  \item{dis}{Dissimilarity data inheriting from class \code{dist} or
    a an object, such as a matrix, that can be converted to a
    dissimilarity matrix. Functions \code{\link{vegdist}} and
    \code{\link{dist}} are some functions producing suitable
    dissimilarity data.}
  \item{toolong}{ Shortest dissimilarity regarded as \code{NA}.
    The function uses a fuzz factor, so
    that dissimilarities close to the limit will be made \code{NA}, too.
    If \code{toolong = 0} (or negative), no dissimilarity is regarded
    as too long.
  }
  \item{trace}{Summarize results of \code{distconnected}}
  \item{x}{Community data.}
  
}
\details{
  Data sets are disconnected if they have sample plots or groups of
  sample plots which share no species with other sites or groups of
  sites. Such data sets cannot be sensibly ordinated by any
  unconstrained method because these subsets cannot be related to each
  other. For instance, correspondence analysis will polarize these
  subsets with eigenvalue 1. Neither can such dissimilarities be
  transformed with \code{\link{stepacross}}, because there is no path
  between all points, and result will contain \code{NA}s. Function
  \code{distconnected} will find such subsets in dissimilarity
  matrices. The function will return a grouping vector that can be used
  for sub-setting the data. If data are connected, the result vector will
  be all \eqn{1}s. The connectedness between two points can be defined
  either by a threshold \code{toolong} or using input dissimilarities
  with \code{NA}s.

  Function \code{no.shared} returns a \code{dist} structure having value
  \code{TRUE} when two sites have nothing in common, and value
  \code{FALSE} when they have at least one shared species. This is a
  minimal structure that can be analysed with \code{distconnected}. The
  function can be used to select dissimilarities with no shared species
  in indices which do not have a fixed upper limit.
  
  Function \code{distconnected} uses depth-first search
  (Sedgewick 1990). 
}
\value{
  Function \code{distconnected} returns a vector for
  observations using integers to identify connected groups. If the data
  are connected, values will be all \code{1}. Function \code{no.shared}
  returns an object of class \code{\link{dist}}.
}
\references{
 Sedgewick, R. (1990). \emph{Algorithms in C}. Addison Wesley. 
}
\author{ Jari Oksanen }

\seealso{\code{\link{vegdist}} or \code{\link{dist}} for getting
  dissimilarities, \code{\link{stepacross}} for a case where you may need
  \code{distconnected}, and for connecting points \code{\link{spantree}}.
}

\examples{
## There are no disconnected data in vegan, and the following uses an
## extremely low threshold limit for connectedness. This is for
## illustration only, and not a recommended practice.
data(dune)
dis <- vegdist(dune)
gr <- distconnected(dis, toolong=0.4)
# Make sites with no shared species as NA in Manhattan dissimilarities
dis <- vegdist(dune, "manhattan")
is.na(dis) <- no.shared(dune)
}
\keyword{ multivariate}

