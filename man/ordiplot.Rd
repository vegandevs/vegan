\name{ordiplot}
\alias{ordiplot}
\alias{identify.ordiplot}
\alias{scores.ordiplot}
\alias{points.ordiplot}
\alias{text.ordiplot}

\title{ Alternative plot and identify Functions for Ordination }

\description{
  Function \code{ordiplot} is an alternative plotting function which
  can be worked with any \pkg{vegan} ordination result and many
  non-\pkg{vegan} results. In addition, \code{plot} functions for
  \pkg{vegan} ordinations return invisibly an \code{"ordiplot"} result
  object, and this allows using \code{ordiplot} support functions with
  this result: \code{identify} can be used to add labels to selected
  site, species or constraint points, and \code{points} and
  \code{text} can add elements to the plot.
}

\usage{
ordiplot(ord, choices = c(1, 2), type="points", display, xlim, ylim,
     cex = 0.7, ...)
\method{identify}{ordiplot}(x, what, labels,  ...)
\method{points}{ordiplot}(x, what, select, arrows = FALSE, ...)
\method{text}{ordiplot}(x, what, labels, select, arrows = FALSE,
    length = 0.05, ...)
}

\arguments{
  \item{ord}{A result from an ordination.}
  \item{choices}{Axes shown. }
  \item{type}{The type of graph which may be \code{"points"},
    \code{"text"} or \code{"none"} for any ordination method.}
  \item{display}{Display only "sites" or "species". The default for most
    methods is to display both, but for \code{\link{cca}},
    \code{\link{rda}} and \code{\link{capscale}} it is the same as in
    \code{\link{plot.cca}}.  }
  \item{xlim, ylim}{the x and y limits (min,max) of the plot.}
  \item{cex}{Character expansion factor for points and text.}
  \item{\dots}{Other graphical parameters. }
  \item{x}{A result object from \code{ordiplot}.}
  \item{what}{Items identified in the ordination plot. The types depend
    on the kind of plot used. Most methods know \code{sites} and
    \code{species}, functions \code{\link{cca}} and \code{\link{rda}}
    know in addition 
    \code{constraints} (for LC scores), \code{centroids} and
    \code{biplot}, and \code{\link{plot.procrustes}} ordination plot has
    \code{heads} and \code{points}.}
  \item{labels}{Optional text used for labels. Row names will be used if
    this is missing.}
  \item{arrows}{Draw arrows from the origin. This will always be
    \code{TRUE} for biplot scores and its value will be ignored. Setting
    this \code{TRUE} will draw arrows for any type of scores. This
    allows, e.g, using biplot arrows for species. The arrow head will be
    at the value of scores, and possible text is moved outwards.}
  \item{length}{Length of arrow heads (see \code{\link{arrows}}).}
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of indices
    of displayed items.}

}
\details{
  Function \code{ordiplot} draws an ordination diagram using black circles for
  sites and red crosses for species.  It returns invisibly an object of
  class \code{ordiplot} which can be used by \code{identify.ordiplot}
  to label selected sites or species, or constraints in
  \code{\link{cca}} and \code{\link{rda}}.

  The function can handle output from several alternative ordination
  methods. For \code{\link{cca}}, \code{\link{rda}} and
  \code{\link{decorana}} it uses their \code{plot} method with option
  \code{type = "points"}. In addition, the \code{plot} functions of
  these methods return invisibly an \code{ordiplot} object which can
  be used by \code{identify.ordiplot} to label points. For other
  ordinations it relies on \code{\link{scores}} to extract the scores.

  For full user control of plots, it is best to call \code{ordiplot}
  with \code{type = "none"} and save the result, and then add sites and
  species using \code{points.ordiplot} or \code{text.ordiplot} which
  both pass all their arguments to the corresponding default graphical
  functions. The functions can be chained with pipes which allows an
  alternative intuitive way of building up plots. 
}
\value{
  
  Function \code{ordiplot} returns invisibly an object of class
  \code{ordiplot} with used scores. In general, \pkg{vegan} \code{plot}
  functions for ordination results will also return an invisible
  \code{ordiplot} object. If the \code{plot(..., type = "n")} was used
  originally, the plot is empty, and items can be added with the
  invisible object. Functions \code{points} and \code{text} will return
  their input object without modification, which allows chaining these
  commands with pipes. Function \code{identify.ordiplot} uses this
  object to label the point.

}

\author{
  Jari Oksanen
}
\note{
  The purpose of these functions is to provide similar functionality as
  the \code{plot}, \code{plotid} and \code{specid} methods in library
  \code{labdsv}. The functions are somewhat limited in parametrization,
  but you can call directly the standard \code{\link{identify}} and
  \code{\link{plot}} functions for a better user control.
}

\seealso{ \code{\link{identify}} for basic operations, \code{\link{plot.cca}},
  \code{\link{plot.decorana}}, \code{\link{plot.procrustes}} which also
  produce objects for
  \code{identify.ordiplot} and \code{\link{scores}} for extracting
  scores from non-\code{vegan} ordinations. 
}

\examples{
## Draw a plot for a non-vegan ordination (cmdscale).
data(dune)
dune.dis <- vegdist(wisconsin(dune))
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, dune, expand = TRUE)
pl <- ordiplot(dune.mds, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "species", col="blue", cex=0.9)
\dontrun{
## same plot using pipes (pipes |> are available from R version 4.1.0)
if (getRversion() >= "4.1") {
ordiplot(dune.mds, type="n") |>
  points("sites", pch=21, col="red", bg="yellow") |>
  text("species", col="blue", cex=0.9)
## Some people think that species should be shown with arrows in PCA.
## Other ordination methods also return an invisible ordiplot object and
## we can use pipes to draw those arrows.
mod <- rda(dune)
plot(mod, type="n") |>
  points("sites", pch=16, col="red") |>
  text("species", arrows = TRUE, length=0.05, col="blue")
}
}
## Default plot of the previous using identify to label selected points
\dontrun{
pl <- ordiplot(dune.mds)
identify(pl, "spec")}
}
\keyword{ hplot }
\keyword{ iplot }
\keyword{ aplot }
