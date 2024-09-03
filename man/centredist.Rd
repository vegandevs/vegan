\name{centredist}
\alias{centredist}
\alias{centredist.betadisper}
\alias{centredist.default}

\title{
  Distances of Points to Class Centroids
}
\description{
  Function finds distances of points to class centroids.
}
\usage{
\method{centredist}{default}(x, centres, distance = c("euclidean",
  "mahalanobis"),
    ...)
\method{centredist}{betadisper}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Matrix for the default method or a\code{betadisper} result object.
  }
  \item{centres}{Factor to define classes for which centres are found.}
  \item{distance}{Use either simple Euclidean distance of Mahalanobis
    distances that take into account the shape of covariance ellipses.}
  \item{\dots}{Other arguments to functions (ignored).}
}
\details{
  See https://github.com/vegandevs/vegan/issues/606
}
\value{
  Function returns a list with components
  \item{centre}{Centre to which a point was originally allocated.}
  \item{nearest}{Nearest centre to a point.}
  \item{distances}{Distance matrix where each observation point is a
    row, and column gives the distance of the point to that centre.}

}

\author{
  Jari Oksanen
}
\note{
  Persons at Work.
}


\seealso{
  \code{\link{betadisper}}.
}
\examples{
data(dune, dune.env)
d <- vegdist(dune)
mod <- with(dune.env, betadisper(d, Management))
centredist(mod)
}

