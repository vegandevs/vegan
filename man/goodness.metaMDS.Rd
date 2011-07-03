\name{goodness.metaMDS}
\alias{goodness.metaMDS}
\alias{goodness.monoMDS}
\alias{stressplot}
\alias{stressplot.default}
\alias{stressplot.monoMDS}

\title{Goodness of Fit and Shepard Plot for Nonmetric Multidimensional Scaling }
\description{
  Function \code{goodness.metaMDS} find goodness of fit measure for
  points in nonmetric multidimensional scaling, and function
  \code{stressplot} makes a \code{\link[MASS]{Shepard}} diagram.
}
\usage{
\method{goodness}{metaMDS}(object, dis, ...)
\method{stressplot}{default}(object, dis, pch, p.col = "blue", l.col = "red", 
    lwd = 2, ...) 
}

\arguments{
  \item{object}{A result object from \code{\link{metaMDS}}, 
    \code{\link{monoMDS}} or \code{\link[MASS]{isoMDS}}. }
  \item{dis}{Dissimilarities.  This should  not used with
    \code{\link{metaMDS}} or \code{\link{monoMDS}}, but must be used with
    \code{\link[MASS]{isoMDS}}.}
  \item{pch}{Plotting character for points.  Default is dependent on the
    number of points. }
  \item{p.col, l.col}{Point and line colours.}
  \item{lwd}{Line width. For \code{\link{monoMDS}} the default is
    \code{lwd = 1} if more than two lines are drawn, and \code{lwd =
    2} otherwise.}
  \item{\dots}{Other parameters to functions, e.g. graphical parameters.}
}
\details{
  Function \code{goodness.metaMDS} finds a goodness of fit statistic for
  observations (points).  This is defined so that sum of squared values
  is equal to squared stress.  Large values indicate poor fit. 

  Function \code{stressplot} draws a Shepard diagram which is a plot
  of ordination distances and monotone or linear fit line against
  original dissimilarities.  In addition, it displays two
  correlation-like statistics on the goodness of fit in the graph.
  The nonmetric fit is based on stress \eqn{S} and defined as \eqn{R^2
  = 1-S^2}{R2 = 1-S*S}.  The \dQuote{linear fit} is the squared
  correlation between fitted values and ordination distances. For
  \code{\link{monoMDS}}, the \dQuote{linear fit} and \eqn{R^2}{R2}
  from \dQuote{stress type 2} are equal.

  Both functions can be used with \code{\link{metaMDS}},
  \code{\link{monoMDS}} and \code{\link[MASS]{isoMDS}}.  The original
  dissimilarities should not be given for \code{\link{monoMDS}} or
  \code{\link{metaMDS}} results (the latter tries to reconstruct the
  dissimilarities using \code{\link{metaMDSredist}} if
  \code{\link[MASS]{isoMDS}} was used as its engine).  With
  \code{\link[MASS]{isoMDS}} the dissimilarities must be given.  In
  either case, the functions inspect that dissimilarities are
  consistent with current ordination, and refuse to analyse
  inconsistent dissimilarities.  Function \code{goodness.metaMDS} is
  generic in \pkg{vegan}, but you must spell its name completely with
  \code{\link[MASS]{isoMDS}} which has no class.
} 

\value{ Function \code{goodness} returns a vector of values. Function
  \code{stressplot} returns invisibly an object with itmes for
  original dissimilarities, ordination distances and fitted values.  }

\author{Jari Oksanen. }

\seealso{\code{\link{metaMDS}},  \code{\link{monoMDS}}, 
  \code{\link[MASS]{isoMDS}}, \code{\link[MASS]{Shepard}}. }

\examples{
data(varespec)
mod <- metaMDS(varespec)
stressplot(mod)
gof <- goodness(mod)
gof
plot(mod, display = "sites", type = "n")
points(mod, display = "sites", cex = 2*gof/mean(gof))
}
\keyword{ multivariate }

