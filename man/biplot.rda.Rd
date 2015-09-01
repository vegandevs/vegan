\name{biplot.rda}
\alias{biplot.rda}

\title{PCA biplot}
\description{
  Draws a PCA biplot with species scores indicated by biplot arrows
}
\usage{
\method{biplot}{rda}(x, choices = c(1, 2), scaling = 2,
       display = c("sites", "species"), type, xlim, ylim, col = c(1,2), 
       const, ...)
}

\arguments{
  \item{x}{A \code{\link{rda}} result object.}
  \item{choices}{Axes to show.}
   \item{scaling}{Scaling for species and site scores. Either species
    (\code{2}) or site (\code{1}) scores are scaled by eigenvalues, and
    the other set of scores is left unscaled, or with \code{3} both are
    scaled symmetrically by square root of eigenvalues. With negative
    scaling values in \code{rda}, species scores are divided by standard
    deviation of each species and multiplied with an equalizing
    constant. Unscaled raw scores stored in the result can be accessed
    with \code{scaling = 0}.}
  \item{display}{Scores shown.  These must some of the alternatives
    \code{"species"} for species scores, and/or \code{"sites"} for site
    scores.}
  \item{type}{Type of plot: partial match to \code{text}
    for text labels, \code{points} for points, and \code{none} for
    setting frames only. If omitted, \code{text} is selected for
    smaller data sets, and \code{points} for larger. Can be of length 2
    (e.g. \code{type = c("text", "points")}), in which case the first
    element describes how species scores are handled, and the second how
    site scores are drawn.}
  \item{xlim, ylim}{the x and y limits (min, max) of the plot.}
  \item{col}{Colours used for sites and species (in this order). If only
    one colour is given, it is used for both.}
  \item{const}{General scaling constant for \code{\link{scores.rda}}.}
  \item{\dots}{Other parameters for plotting functions.}
}
\details{
  Produces a plot or biplot of the results of a call to
  \code{\link{rda}}. It is common for the "species" scores in a PCA to
  be drawn as biplot arrows that point in the direction of increasing
  values for that variable. The \code{biplot.rda} function provides a
  wrapper to \code{plot.cca} to allow the easy production of such a
  plot.

  \code{\link{biplot.rda}} is only suitable for unconstrained models. If
  used on an ordination object with constraints, an error is issued.

  If species scores are drawn using \code{"text"}, the arrows are drawn
  from the origin to 0.85 * species score, whilst the labels are
  drawn at the species score. If the type used is \code{"points"}, then
  no labels are drawn and therefore the arrows are drawn from the origin
  to the actual species score.
}
\value{
  The \code{plot} function returns invisibly a plotting structure which
  can be used by \code{\link{identify.ordiplot}} to identify
  the points or other functions in the \code{\link{ordiplot}} family.
}
\author{Gavin Simpson, based on \code{\link{plot.cca}} by Jari Oksanen.}
\seealso{\code{\link{plot.cca}}, \code{\link{rda}} for something to
  plot, \code{\link{ordiplot}} for an alternative plotting routine
  and more support functions, and \code{\link{text}},
  \code{\link{points}} and \code{\link{arrows}} for the basic routines.}
\examples{
data(dune)
mod <- rda(dune, scale = TRUE)
biplot(mod, scaling = 3)

## different type for species and site scores
biplot(mod, scaling = 3, type = c("text", "points"))
}
\keyword{hplot}
