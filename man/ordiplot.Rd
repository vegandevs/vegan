\name{ordiplot}
\alias{ordiplot}
\alias{identify.ordiplot}
\alias{scores.ordiplot}
\alias{points.ordiplot}
\alias{text.ordiplot}

\title{ Alternative plot and identify Functions for Ordination }

\description{ Function \code{ordiplot} is an alternative plotting
  function which works with any \pkg{vegan} ordination object and many
  non-\pkg{vegan} objects. In addition, \code{plot} functions for
  \pkg{vegan} ordinations return invisibly an \code{"ordiplot"} object,
  and this allows using \code{ordiplot} support functions with this
  result: \code{identify} can be used to add labels to selected site,
  species or constraint points, and \code{points} and \code{text} can
  add elements to the plot, and  used in a pipe to add scores into plot
  by layers.  }

\usage{
ordiplot(ord, choices = c(1, 2), type="points", display, optimize = FALSE,
    arrows = FALSE, length = 0.05, arr.mul, xlim, ylim, ...)
\method{points}{ordiplot}(x, what, select, arrows = FALSE,
    length = 0.05, arr.mul, ...)
\method{text}{ordiplot}(x, what, labels, select, optimize = FALSE,
    arrows = FALSE, length = 0.05, arr.mul, bg, ...)
\method{identify}{ordiplot}(x, what, labels,  ...)
}

\arguments{
  \item{ord}{A result from an ordination.}
  \item{choices}{Axes shown. }
  \item{type}{The type of graph which may be \code{"points"},
    \code{"text"} or \code{"none"} for any ordination method.}
  \item{display}{Display only "sites" or "species". The default for most
    methods is to display both, but for \code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}} and \code{\link{capscale}}
    it is the same as in \code{\link{plot.cca}}.  }
  \item{xlim, ylim}{the x and y limits (min,max) of the plot.}
  \item{\dots}{Other graphical parameters. }
  \item{x}{A result object from \code{ordiplot}.}
  \item{what}{Items identified in the ordination plot. The types depend
    on the kind of plot used. Most methods know \code{sites} and
    \code{species}, functions \code{\link{cca}} and \code{\link{rda}}
    know in addition 
    \code{constraints} (for LC scores), \code{centroids}, \code{biplot}
    and \code{regression}, and \code{\link{plot.procrustes}} ordination
    plot has \code{heads} and \code{points}.}
  \item{labels}{Optional text used for labels. Row names of scores will
    be used if this is missing. If \code{select} is used, \code{labels} 
    are given only to selected items in the order they occur in the scores.}
  \item{optimize}{Optimize locations of text to reduce overlap and plot
    point in the actual locations of the scores. Uses
    \code{\link{ordipointlabel}}. }
   \item{arrows}{Draw arrows from the origin. This will always be
    \code{TRUE} for biplot and regression scores in constrained
    ordination (\code{\link{cca}} etc.). Setting this \code{TRUE} will
    draw arrows for any type of scores. This allows, e.g, using biplot
    arrows for species. The arrow head will be at the value of scores,
    and possible text is moved outwards.}
  \item{length}{Length of arrow heads (see \code{\link{arrows}}).}
  \item{arr.mul}{Numeric multiplier to arrow lenghts; this will also set
    \code{arrows = TRUE}. The default is to automatically adjust arrow
    lengths with \code{"biplot"} and \code{"regression"} scores and else
    use unmodified scores.}
  \item{bg}{Background colour for labels. If \code{bg} is set, the
    labels are displayed with \code{\link{ordilabel}} instead of
    \code{\link{text}}.  }
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of indices
    of displayed items.}

}
\details{
  
  Function \code{ordiplot} draws an ordination diagram with default of
  black circles for sites and red crosses for species.  It returns
  invisibly an object of class \code{ordiplot}.

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
  functions. Alternatively, \code{points} and \code{text} can be used in
  pipe which allows an intuitive way of building up plots by layers. In
  addition, function \code{\link{ordilabel}} and
  \code{\link{ordipointlabel}} can be used in pipe after \code{ordiplot}
  or other \pkg{vegan} ordination \code{plot} commands. See Examples. 
}
\value{
  
  Function \code{ordiplot} returns invisibly an object of class
  \code{ordiplot} with used scores. In general, \pkg{vegan} \code{plot}
  functions for ordination results will also return an invisible
  \code{ordiplot} object. If the \code{plot(..., type = "n")} was used
  originally, the plot is empty, and items can be added with the
  invisible object. Functions \code{points} and \code{text} return
  their input object without modification, which allows chaining these
  commands with pipes. Function \code{identify.ordiplot} uses this
  object to label the point.

}

\author{
  Jari Oksanen
}

\seealso{
  With argument \code{bg} function calls \code{\link{ordilabel}}
  to draw text on non-transparent label, and with argument
  \code{optimize = TRUE} function calls \code{\link{ordipointlabel}} to
  optimize the locations of text labels to minimize
  over-plotting. Functions \code{\link{ordilabel}} and
  \code{\link{ordipointlabel}} can be used in a pipe together with
  \code{ordiplot} methods \code{text} and \code{points}. Function
  \code{\link{plot.cca}} uses \code{ordiplot} methods \code{text} and
  \code{points} in configurable plots, and these accept the
  arguments of the \code{ordiplot} methods described here.
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

## same plot using pipes (|>)
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

## Default plot of the previous using identify to label selected points
\dontrun{
pl <- ordiplot(dune.mds)
identify(pl, "spec")}
}
\keyword{ hplot }
\keyword{ iplot }
\keyword{ aplot }
