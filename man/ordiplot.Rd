\name{ordiplot}
\alias{ordiplot}
\alias{identify.ordiplot}
\alias{scores.ordiplot}
\alias{points.ordiplot}
\alias{text.ordiplot}

\title{ Alternative plot and identify Functions for Ordination }
\description{
  Ordination plot function especially for congested plots. Function
  \code{ordiplot} always plots only unlabelled points, but
  \code{identify.ordiplot} can be used to add labels to selected site,
  species or constraint points.  Function \code{identify.ordiplot} can
  be used to identify points from \code{\link{plot.cca}},
  \code{\link{plot.decorana}}, \code{\link{plot.procrustes}} and
  \code{\link{plot.rad}} as well.
}
\usage{
ordiplot(ord, choices = c(1, 2), type="points", display, xlim, ylim,
     cex = 0.7, ...)
\method{identify}{ordiplot}(x, what, labels,  ...)
\method{points}{ordiplot}(x, what, select, ...)
\method{text}{ordiplot}(x, what, labels, select, ...)
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
    \code{constraints} (for `LC' scores), \code{centroids} and
    \code{biplot}, and \code{\link{plot.procrustes}} ordination plot has
    \code{heads} and \code{points}.}
  \item{labels}{Optional text used for labels. Row names will be used if
  this is missing.}
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
  functions.
}
\value{
  Function \code{ordiplot} returns invisibly an object of class
  \code{ordiplot} with items \code{sites}, \code{species} and
  \code{constraints} (if these are available in the ordination
  object). Function \code{identify.ordiplot} uses this object to label
  the point.
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
# Draw a plot for a non-vegan ordination (cmdscale).
data(dune)
dune.dis <- vegdist(wisconsin(dune))
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, dune, expand = TRUE)
fig <- ordiplot(dune.mds, type = "none")
points(fig, "sites", pch=21, col="red", bg="yellow")
text(fig, "species", col="blue", cex=0.9)
# Default plot of the previous using identify to label selected points
\dontrun{
fig <- ordiplot(dune.mds)
identify(fig, "spec")}

}
\keyword{ hplot }
\keyword{ iplot }
\keyword{ aplot }
