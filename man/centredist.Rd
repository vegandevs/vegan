\name{centredist}
\alias{centredist}
\alias{centredist.betadisper}

\title{
  Distances of Points to Class Centroids
}
\description{
  Function finds distances of points to class centroids.
}
\usage{
\method{centredist}{betadisper}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    \code{betadisper} result object
}
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

