\name{goodness.metaMDS}
\alias{goodness.metaMDS}
\alias{stressplot}

\title{Goodness of Fit and Shepard Plot for Nonmetric Multidimensional Scaling }
\description{
  Function \code{goodness.metaMDS} find goodness of fit measure for
  points in nonmetric multidimensional scaling, and function
  \code{stressplot} makes a \code{\link[MASS]{Shepard}} diagram.
}
\usage{
\method{goodness}{metaMDS}(object, dis, ...)
stressplot(object, dis, pch, p.col = "blue", l.col = "red", lwd = 2, 
    ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A result object from \code{\link{metaMDS}} or
    \code{\link[MASS]{isoMDS}}. }
  \item{dis}{Dissimilarities.  Normally this should  not used with
    \code{\link{metaMDS}}, but should be always used with
    \code{\link[MASS]{isoMDS}}.}
  \item{pch}{Plotting character for points.  Default is dependent on the
    number of points. }
  \item{p.col, l.col}{Point and line colours.}
  \item{lwd}{Line width.}
  \item{\dots}{Other parameters to functions, e.g. graphical parameters. }
}
\details{
  Function \code{goodness.metaMDS} finds a goodness of fit statistic for
  observations (points).  This is defined so that sum of squared values
  is equal to squared stress.  Large values indicate poor fit.

  Function \code{stressplot} is a wrapper to \code{\link[MASS]{Shepard}}
  function in \pkg{MASS} package.  It plots ordination distances
  against original dissimilarities, and draws a step line of the
  nonlinear fit.  In addition, it adds to the graph two correlation-like
  statistics on the goodness of fit.  The nonmetric fit is based on stress \eqn{S} and
  defined as \eqn{\sqrt{1-S^2}}{sqrt(1-S^2)}.  The \dQuote{linear fit}
  is the correlation between fitted values and ordination distances.

  Both functions can be used both with \code{\link{metaMDS}} and with
  \code{\link[MASS]{isoMDS}}.  With \code{\link{metaMDS}}, the functions
  try to reconstruct the dissimilarities using
  \code{\link{metaMDSredist}}, and dissimilarities should not be given.
  With  \code{\link[MASS]{isoMDS}} the dissimilarities must be given.
  In either case, the functions inspect that dissimilarities are
  consistent with current ordination, and refuse to analyse inconsistent
  dissimilarities.  Function \code{goodness.metaMDS} is generic in
  \code{vegan}, but you must spell its name completely with
  \code{\link[MASS]{isoMDS}} which has no class.
}
\value{
  Function \code{goodness} returns a vector of values. Function
  \code{stressplot} returns invisibly a \code{\link[MASS]{Shepard}}
  object. 
}
\author{Jari Oksanen. }

\seealso{\code{\link{metaMDS}},  \code{\link[MASS]{isoMDS}},
  \code{\link[MASS]{Shepard}}. }
\examples{
data(varespec)
mod <- metaMDS(varespec)
stressplot(mod)
gof <- goodness(mod)
gof
plot(mod, display = "sites", type = "n")
points(mod, display = "sites", cex = gof/2)
}
\keyword{ multivariate }

