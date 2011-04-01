\name{monoMDS}
\alias{monoMDS}

\title{
  Kruskal's Non-metric Multidimensional Scaling
}
\description{
  Function implements Kruskal's (1964a,b) non-metric multidimensional scaling 
  using monotone regression and primary (\dQuote{weak}) treatment of ties.
}
\usage{
monoMDS(dist, y, k = 2, model = c("global", "local", "hybrid"), 
    threshold = 0.8, maxit = 200, tol = 1e-04, ...)
}

\arguments{
  \item{dist}{Input dissimilarities.}
  \item{y}{Starting configuration. A random configuration will be
    generated if this is missing.}
  \item{k}{Number of dimensions.}
  \item{model}{MDS model: \code{"global"} is normal non-metric MDS
    with a monotone regression, \code{"local"} is non-metric MDS with
    separate regressions for each point, and \code{"hybrid"} uses
    linear MDS for dissimilarities below the threshold, and global
    non-metric MDS for dissimilarities above the threshold. }
  \item{threshold}{Dissimilarity above which monotone regression is
    used instead of linear regression. }
  \item{maxit}{Maximum number of iterations.}
  \item{tol}{Convergence tolerance.}
  \item{\dots}{Other parameters to the functions (ignored).}
}

\details{
  This is the first version of the function, and still under construction.
}

\value{
  A complicated object (under construction).
}

\references{
  Kruskal, J.B. 1964a. Multidimensional scaling by optimizing
      goodness-of-fit to a nonmetric hypothesis. \emph{Psychometrika}
      29, 1--28.

  Kruskal, J.B. 1964b. Nonmetric multidimensional scaling: a numerical
      method. \emph{Psychometrika} 29, 115--129.
}

\author{
Peter R. Michin (Fortran core) and Jari Oksanen (R interface).
}

\note{
   Still under development.
}

\seealso{ \code{\link[vegan]{metaMDS}}, \code{\link[MASS]{isoMDS}},
  \code{\link[smacof]{smacofSym}}, \code{\link[ecodist]{nmds}}. }

\examples{
##
}
\keyword{ multivariate }

