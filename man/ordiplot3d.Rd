\name{ordiplot3d}
\alias{ordiplot3d}
\alias{ordirgl}
\alias{orglpoints}
\alias{orgltext}
\alias{orglsegments}
\alias{orglspider}

\title{Three-Dimensional and Dynamic Ordination Graphics }
\description{
  Function \code{ordiplot3d} displays three-dimensional ordination
  graphics using \code{\link[scatterplot3d]{scatterplot3d}}.  Function
  \code{ordirgl} displays three-dimensional dynamic ordination graphs
  which can be rotated and zoomed into using \code{\link[rgl]{rgl}}
  package. Both work with all ordination
  results form \code{vegan} and all ordination results known by
  \code{\link{scores}} function. 
}
\usage{
ordiplot3d(object, display = "sites", choices = 1:3, ax.col = 2,
        arr.len = 0.1, arr.col = 4, envfit, xlab, ylab, zlab, ...)
ordirgl(object, display = "sites", choices = 1:3, type = "p", 
        ax.col = "red", arr.col = "yellow", text, envfit, ...)
orglpoints(object, display = "sites", choices = 1:3, ...)
orgltext(object, text, display = "sites", choices = 1:3, justify = "center", 
        adj = 0.5, ...)
orglsegments(object, groups, display = "sites", choices = 1:3, ...)
orglspider(object, groups, display = "sites", w = weights(object, display),
        choices = 1:3, ...)
}
\arguments{
  \item{object}{An ordination result or any object known by \code{\link{scores}}. }
  \item{display}{Display \code{"sites"} or \code{"species"} or other
    ordination object recognized by \code{\link{scores}}. }
  \item{choices}{Selected three axes. }
    \item{arr.len}{'Length' (width) of arrow head passed to
    \code{\link{arrows}} function. }
  \item{arr.col}{Colour of biplot \code{\link{arrows}} and centroids of
    environmental variables. }
  \item{type}{The type of plots: \code{"p"} for points or \code{"t"} for
    text labels.}
  \item{ax.col}{Axis colour (concerns only the crossed axes through the
    origin).}
  \item{text}{Text to override the default with \code{type =
    "t"}.}
  \item{envfit}{Fitted environmental variables from \code{\link{envfit}}
    displayed in the graph.}
  \item{xlab, ylab, zlab}{Axis labels passed to
    \code{\link[scatterplot3d]{scatterplot3d}}. If missing, labels  are
    taken from the
    ordination result. Set to \code{NA} to suppress labels. }
  \item{justify, adj}{Text justification passed to
    \code{\link[rgl]{rgl.texts}}. One of these is used depending on the
    version of \pkg{rgl} installed.}
  \item{groups}{Factor giving the groups for which the graphical item is
    drawn.}
  \item{w}{Weights used to find the average within group. Weights are
    used automatically for \code{\link{cca}}
    and \code{\link{decorana}} results, unless undone by the
    user. \code{w=NULL} sets equal weights to all points. }
  \item{\dots}{Other parameters passed to graphical functions. }
}
\details{
  Both function display three-dimensional ordination graphics. Function
  \code{ordiplot3d} plots static scatter diagrams using
  \code{\link[scatterplot3d]{scatterplot3d}}. Function \code{ordirgl}
  plots dynamic graphics using OpenGL  in \code{\link[rgl]{rgl}}. Both
  functions use most default settings of underlying graphical functions,
  and you must consult their help pages to change graphics to suit your
  taste (see \code{\link[scatterplot3d]{scatterplot3d}},
  \code{\link[rgl]{rgl}},
  \code{\link[rgl]{rgl.points}},\code{\link[rgl]{rgl.texts}}). Both
  functions will display only one selected set of \code{\link{scores}},
  typically either \code{"sites"} or \code{"species"}, but
  for instance \code{\link{cca}} also has \code{"lc"} scores. In
  constrained ordination (\code{\link{cca}}, \code{\link{rda}},
  \code{\link{capscale}}), biplot arrows and centroids are always
  displayed similarly as in two-dimensional plotting function
  \code{\link{plot.cca}}.  Alternatively, it is possible to display
  fitted environmental vectors or class centroids from
  \code{\link{envfit}} in both graphs.  These are displayed similarly as
  the results of constrained ordination, and they can be shown only for
  non-constrained ordination. The user must remember to specify at least
  three axes in \code{\link{envfit}} if the results are used with these
  functions. 

  Function \code{ordiplot3d} plots only points. However, it returns
  invisibly an object inheriting from \code{\link{ordiplot}} so that you
  can use \code{\link{identify.ordiplot}} to identify \code{"points"} or
  \code{"arrows"}. The underlying
  \code{\link[scatterplot3d]{scatterplot3d}} function accepts \code{type
    = "n"} so that only the axes, biplot arrows and centroids of
  environmental variables will be plotted, and the ordination scores can
  be added with \code{\link{text.ordiplot}} or
  \code{\link{points.ordiplot}}. Further, you can use any functions from
  the \code{\link{ordihull}} family with the invisible result of
  \code{\link{ordiplot3d}}, but you must remember to specify the
  \code{display} as \code{"points"} or \code{"arrows"}. To change the
  viewing angle, orientation etc. you must see
  \code{\link[scatterplot3d]{scatterplot3d}}.

  Function \code{ordigl} makes a dynamic three-dimensional graph that
  can be rotated with mouse, and zoomed into with mouse buttons or wheel
  (but Mac users with one-button mouse should see
  \code{\link[rgl]{rgl.viewpoint}}), or try ctrl-button. MacOS X users
  must start \code{X11} before calling \code{\link[rgl]{rgl}} commands.
  Function \code{ordirgl} uses default settings, and you should consult the
  underlying functions  \code{\link[rgl]{rgl.points}},
  \code{\link[rgl]{rgl.texts}} to see how to control the
  graphics. Function \code{ordirgl} always cleans its graphic window
  before drawing.  Functions \code{orglpoints} adds points and
  \code{orgltext} adds text to existing \code{ordirgl} windows.  In
  addition, function \code{orglsegments} combines points within
  \code{"groups"} with line segments similarly as
  \code{\link{ordisegments}}. Function \code{orglspider} works similarly
  as \code{\link{ordispider}}: it connects points to their weighted
  centroid within \code{"groups"}, and in constrained ordination it can
  connect \code{"wa"} or weighted averages scores to corresponding
  \code{"lc"} or linear combination scores if \code{"groups"} is
  missing. In addition, basic \code{rgl} functions 
  \code{\link[rgl]{rgl.points}}, \code{\link[rgl]{rgl.texts}},
 \code{\link[rgl]{rgl.lines}} and many others can be used. 
}
 \value{
  Function \code{ordiplot3d} returns invisibly an object of class
  \code{"ordiplot3d"} inheriting from \code{\link{ordiplot}}. The return
  object will contain the coordinates projected onto two dimensions for
  \code{"points"}, and possibly for the heads of \code{"arrows"} and
  \code{"centroids"} of environmental variables. Functions like
  \code{\link{identify.ordiplot}}, \code{\link{points.ordiplot}},
  \code{\link{text.ordiplot}} can use this result, as well as
  \code{\link{ordihull}} and other functions documented with the
  latter. In addition, the result will contain the object returned by
  \code{\link[scatterplot3d]{scatterplot3d}}, including function
  \code{xyz.converter} which projects three-dimensional
  coordinates onto the plane used in the current plot. 
  Function \code{ordirgl} returns nothing.
}

\author{Jari Oksanen }

\section{Warning}{Function \code{ordirgl} uses OpenGL package
  \code{\link[rgl]{rgl}}
  which may not be functional in all platforms, and can crash R in some:
  use \code{\link{save.image}} before trying \code{ordirgl}.
  Mac users must start \code{X11} (and first install \code{X11} and some other
  libraries) before being able to use \code{\link[rgl]{rgl}}. It seems
  that \code{\link[rgl]{rgl.texts}} does  not always position the text
  like supposed, and it may be safe to verify text location with
  corresponding points.
}

\note{The user interface of \pkg{rgl} changed in version 0.65, but
  the \code{ordirgl} functions do not yet fully use the new
  capabilities.  However, they should work both in old and new versions
  of \pkg{rgl}. 
  }

\seealso{   \code{\link[scatterplot3d]{scatterplot3d}},
  \code{\link[rgl]{rgl}}, \code{\link[rgl]{rgl.points}},
  \code{\link[rgl]{rgl.texts}}, \code{\link[rgl]{rgl.viewpoint}},
  \code{\link{ordiplot}}, \code{\link{identify.ordiplot}}, 
  \code{\link{text.ordiplot}}, \code{\link{points.ordiplot}},
  \code{\link{ordihull}}, \code{\link{plot.cca}}, \code{\link{envfit}}.
}
\examples{
## Examples are not run, because they need non-standard packages
## 'scatterplot3d' and 'rgl' (and the latter needs user interaction).
#####
#### Default 'ordiplot3d'
\dontrun{
data(dune)
data(dune.env)
ord <- cca(dune ~ A1 + Moisture, dune.env)
ordiplot3d(ord)
#### A boxed 'pin' version
ordiplot3d(ord, type = "h")
#### More user control
pl <- ordiplot3d(ord, angle=15, type="n")
points(pl, "points", pch=16, col="red", cex = 0.7)
#### identify(pl, "arrows", col="blue") would put labels in better positions
text(pl, "arrows", col="blue", pos=3)
text(pl, "centroids", col="blue", pos=1, cex = 1.2)
#### ordirgl
ordirgl(ord, size=2)
ordirgl(ord, display = "species", type = "t")
rgl.quit()
}
}
\keyword{ hplot }
\keyword{ dynamic }
