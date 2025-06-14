\name{biplot.rda}
\alias{biplot.rda}
\alias{biplot.cca}

\title{PCA biplot}
\description{
  Draws a PCA biplot with species scores indicated by biplot arrows
}
\usage{
\method{biplot}{rda}(x, choices = c(1, 2), scaling = "species",
       display = c("sites", "species"), type, xlim, ylim, col = c(1,2), 
       const, correlation = FALSE, ...)
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
    with \code{scaling = 0}.

    The type of scores can also be specified as one of \code{"none"},
    \code{"sites"}, \code{"species"}, or \code{"symmetric"}, which
    correspond to the values \code{0}, \code{1}, \code{2}, and \code{3}
    respectively. Argument \code{correlation} can be used in combination
    with these character descriptions to get the corresponding negative
    value.
  }
  \item{correlation}{logical; if \code{scaling} is a character
    description of the scaling type, \code{correlation} can be used to
    select correlation-like scores for PCA. See argument \code{scaling}
    for details.}
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

  Arrow heads are at the actual species scores, and the possible text is
  projected to direction of the arrow.

}
\value{
  The \code{plot} function returns invisibly a plotting structure which
  can be used by \code{\link{identify.ordiplot}} to identify
  the points or other functions in the \code{\link{ordiplot}} family.
}
\author{Gavin Simpson and Jari Oksanen.}

\note{Prior to \pkg{vegan} 2.7-2 the species \code{text} was at the
  actual location and arrows were shorter than with \code{points}.}

\seealso{\code{\link{plot.cca}} can also draw biplot arrows since
  \pkg{vegan} 2.7-0.}
\examples{
data(dune)
mod <- rda(dune, scale = TRUE)
biplot(mod, scaling = "symmetric")

## plot.cca can do the same
plot(mod, scaling = "symmetric", spe.par = list(arrows=TRUE))

## different type for species and site scores
biplot(mod, scaling = "symmetric", type = c("text", "points"))

## We can use ordiplot pipes to build similar plots with flexible
## control
plot(mod, scaling = "symmetric", type="n") |>
   points("sites", cex=0.7) |>
   text("species", arrows=TRUE, length=0.05, col=2, cex=0.7, font=3)
}

\keyword{hplot}
