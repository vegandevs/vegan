\name{centredist}
\alias{centredist}
\alias{centerdist}
\alias{centredist.betadisper}
\alias{centredist.dbrda}
\alias{centredist.default}
\alias{centredist.cca}
\alias{centredist.rda}
\alias{centredist.wcmdscale}

\title{
  Distances of Points to Class Centroids
}
\description{
  Function finds Euclidean or squared Mahalanobis distances of points
  to class centroids in \code{\link{betadisper}} or ordination results.
}

\usage{
\method{centredist}{default}(x, group, distance = c("euclidean",
    "mahalanobis"), display = "sites", w, ...)
\method{centredist}{betadisper}(x, ...)
\method{centredist}{cca}(x, group, distance = c("euclidean", "mahalanobis"),
    display = c("sites", "species"), rank = 2, ...)
\method{centredist}{rda}(x, group, distance = c("euclidean", "mahalanobis"),
    display = "sites", rank = 2, ...)
\method{centredist}{dbrda}(x, group, distance = c("euclidean","mahalanobis"),
    display = "sites", rank = 2, ...)
\method{centredist}{wcmdscale}(x, group, distance= c("euclidean","mahalanobis"),
    display = "sites", rank = 2, ...) 
}

\arguments{
  \item{x}{An ordination result object, or a matrix for the default
    method or a \code{betadisper} result object.}
  \item{group}{Factor to define classes for which centroids are found.}
  \item{distance}{Use either simple Euclidean distance or squared Mahalanobis
    distances that take into account the shape of covariance ellipses.}
  \item{display}{Kind of \code{\link{scores}}.}
  \item{w}{Weights. If \code{x} has weights (such as \code{\link{cca}}
    result), these will be used.}
  \item{rank}{Number of axes used in ordination methods.}
  \item{\dots}{Other arguments to functions (ignored).}
}

\details{

  Function finds either simple Euclidean distances of all points to
  their centroid as given by argument \code{group} or squared
  Mahalanobis distances of points to the centroid. In addition to the
  distances, it returns information of the original group and the
  nearest group for each point.

  Currently the distances are scaled so that the sum of squared
  Euclidean distance will give the residual (unexplained) eigenvalue
  in constrained ordination with given group. The distances may not
  correspond to the distances in ordination plots with specific
  \code{scaling} and setting argument \code{scaling} gives an
  error. Due to restriction of this scaling, Euclidean distances for
  species scores are only allowed in \code{\link{ca}} and
  \code{\link{cca}}. These restrictions may be relaxed in the future
  versions of the method.

  Squared Mahalanobis distances (\code{\link{mahalanobis}}) are scaled
  by the within group covariance matrix, and they are consistent with
  \code{\link{ordiellipse}}.

  Mahalanobis distances can only be calculated for \code{group} that
  are larger than \code{rank}. For \code{rank = 2}, minimum
  permissible \code{group} has three points. However, even this can
  fail in degenerate cases, such as three points on a line. Euclidean
  distances can be calculated for any \code{rank}. With \code{rank = "full"}
  distance-based methods \code{\link{wcmdscale}}, \code{\link{pco}},
  \code{\link{dbrda}} and \code{\link{betadisper}} will handle
  negative eigenvalues and associated eigenvectors, but with numeric
  rank they will only use real eigenvectors with positive eigenvalues.

  There are specific methods for some ordination models, but many will
  be handled by the default method (including \code{\link{metaMDS}}
  and \code{\link{monoMDS}}).

  The specific method for \code{\link{betadisper}} ignores most
  arguments and calculates the distances as specified in model
  fitting. The \code{\link{betadisper}} result with argument
  \code{type = "centroid"} is equal to \code{\link{wcmdscale}},
  \code{\link{pco}} or \code{\link{dbrda}} with the same
  \code{group} and \code{rank = "full"}. 

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

\note{\code{centerdist} is a synonym of \code{centredist}.}

\examples{
  data(dune, dune.env)
  d <- vegdist(dune) # Bray-Curtis
  mod <- with(dune.env, betadisper(d, Management))
  centredist(mod)
  ## equal to metric scaling when using centroids instead of medians
  modc <- with(dune.env, betadisper(d, Management, type = "centroid"))
  ord <- pco(d)
  all.equal(centredist(modc), centredist(ord, dune.env$Management, rank="full"),
     check.attributes = FALSE)
}


\seealso{
  \code{\link{betadisper}}, \code{\link{goodness.cca}},
  \code{\link{mahalanobis}} and various ordination methods.
}

