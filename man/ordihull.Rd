\name{ordihull}
\alias{ordihull}
\alias{ordiarrows}
\alias{ordisegments}
\alias{ordigrid}
\alias{ordispider}
\alias{ordiellipse}
\alias{ordicluster}
\alias{weights.cca}
\alias{weights.rda}
\alias{weights.decorana}
\alias{summary.ordihull}
\alias{summary.ordiellipse}

\title{Add Graphical Items to Ordination Diagrams}
\description{
  Functions to add convex hulls, arrows, line segments, regular grids of
  points, `spider' graphs, ellipses or cluster dendrogram 
  to ordination diagrams. The
  ordination diagrams can be produced by \code{vegan}
  \code{\link{plot.cca}}, \code{\link{plot.decorana}} or
  \code{\link{ordiplot}}.
}
\usage{
ordihull(ord, groups, display = "sites", draw = c("lines","polygon", "none"),
         show.groups, label = FALSE,  ...)
ordiellipse(ord, groups, display="sites", kind = c("sd","se"), conf,
         draw = c("lines","polygon", "none"), w = weights(ord, display),
         show.groups, label = FALSE, ...)
ordispider(ord, groups, display="sites", w = weights(ord, display),
         show.groups, ...)
ordiarrows(ord, groups, levels, replicates, display = "sites",
         show.groups, startmark, ...)
ordisegments(ord, groups, levels, replicates, display = "sites",
         show.groups, ...)
ordigrid(ord, levels, replicates, display = "sites",  ...)
ordicluster(ord, cluster, prune = 0, display = "sites",
         w = weights(ord, display), ...)
}

\arguments{
  \item{ord}{An ordination object or an \code{\link{ordiplot}} object. }
  \item{groups}{Factor giving the groups for which the graphical item is
    drawn. }
  \item{levels, replicates}{Alternatively, regular
    groups can be defined with arguments \code{levels} and
    \code{replicates}, where \code{levels} gives the number of groups,
    and \code{replicates} the number of successive items at the same
    group.}
  \item{display}{Item to displayed. }
  \item{draw}{Use either \code{\link{lines}} or \code{\link{polygon}} to
    draw the
    line. Graphical parameters are passed to both. The main difference
    is that \code{polygon}s may be filled and non-transparent. With \code{none}
    nothing is drawn, but the function returns the \code{\link{invisible}}
    plotting data.}
  \item{show.groups}{Show only given groups. This can be a vector, or
    \code{TRUE} if you want to show items for which condition is
    \code{TRUE}. This argument makes it possible to use different
    colours and line types for groups. The default is to show all
    groups. }
  \item{label}{Label the centre of the points in hull or the ellipse with 
    the group name.}
  \item{startmark}{plotting character used to mark the first
    item. The default is to use no mark, and for instance, 
    \code{startmark = 1} will draw a circle.  For other plotting characters,
    see \code{pch} in \code{\link{points}}. }
  \item{w}{Weights used to find the average within group. Weights are
    used automatically for \code{\link{cca}}
    and \code{\link{decorana}} results, unless undone by the
    user. \code{w=NULL} sets equal weights to all points. }
  \item{kind}{Whether standard deviations of points (\code{sd}) or
    standard deviations of their (weighted) averages (\code{se}) are
    used.}
  \item{conf}{Confidence limit for ellipses, e.g. 0.95. If given, the
    corresponding \code{sd} or \code{se} is multiplied with the
    corresponding value found from the Chi-squared distribution with
    2df. }
  \item{cluster}{Result of hierarchic cluster analysis, such as
    \code{\link{hclust}} or \code{\link[cluster]{agnes}}.}
  \item{prune}{Number of upper level hierarchies removed from the
    dendrogram. If \code{prune} \eqn{>0}, dendrogram will be
    disconnected.}
   \item{\dots}{Parameters passed to graphical functions such
    as \code{\link{lines}}, \code{\link{segments}},
    \code{\link{arrows}}, \code{\link{polygon}} or to
    \code{\link{scores}} to select axes and
    scaling etc. }
  }
\details{
  Function \code{ordihull} draws \code{\link{lines}} or
  \code{\link{polygon}}s for the convex
  hulls found by function \code{\link{chull}} encircling
  the items in the groups. 
  
  Function \code{ordiellipse} draws \code{\link{lines}} or
  \code{\link{polygon}}s for dispersion \code{\link[ellipse]{ellipse}}
  using either standard deviation of point scores or standard error of
  the (weighted) average of scores, and the (weighted) correlation
  defines the direction of the principal axis of the ellipse. 
  An ellipsoid hull can be drawn with function
  \code{\link[cluster]{ellipsoidhull}} of package \pkg{cluster}.

  Functions \code{ordihull} and \code{ordiellipse} return the
  \code{\link{invisible}} plotting structure. In \code{ordihull} this
  is a list of coordinates of the hull and in \code{ordiellipse} a
  list of covariance matrices and scales used in drawing the
  ellipses. These result objects have a summary method that returns
  the coordinates of the centres of the ellipses or hulls and their
  surface areas in user units. The centres of the hulls may differ
  from the location of the \code{label} which is the centre of the
  points instead of the centre of the polygon. With \code{draw = "none"}
  only the result object is returned and nothing is drawn.

  Function \code{ordiarrows} draws
  \code{\link{arrows}} and \code{ordisegments} draws line
  \code{\link{segments}} between successive items in the
  groups. Function \code{ordigrid} draws line
  \code{\link{segments}} both within the groups and for the
  corresponding items among the groups.

  Function \code{ordispider} draws a `spider' diagram where each point
  is connected to the group centroid with
  \code{\link{segments}}. Weighted centroids are used in the
  correspondence analysis  methods \code{\link{cca}} and
  \code{\link{decorana}} or if the user gives the weights in the
  call. If \code{ordispider} is called with \code{\link{cca}} or
  \code{\link{rda}} result without \code{groups} argument, the function
  connects each `WA' scores to the corresponding `LC' score.

  Function \code{ordicluster} overlays a cluster dendrogram onto
  ordination. It needs the result from a hierarchic clustering such as
  \code{\link{hclust}} or \code{\link[cluster]{agnes}}, or other
  with a similar structure. Function \code{ordicluster} connects
  cluster centroids to each other with line
  \code{\link{segments}}. Function uses centroids of all points in the 
  clusters, and is therefore similar to average linkage methods.
}

\note{These functions add graphical items to ordination graph: You must
  draw a graph first.
  }
\author{ Jari Oksanen }

\seealso{The functions pass parameters to basic graphical functions, and
  you may wish to change the default values in \code{\link{arrows}},
  \code{\link{lines}}, \code{\link{segments}} and
  \code{\link{polygon}}. You can pass
  parameters to \code{\link{scores}} as well. Underlying function for
  \code{ordihull} is \code{\link{chull}}.
}


\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ Management, dune.env)
attach(dune.env)
## pass non-graphical arguments without warnings
plot(mod, type="n", scaling = 3)
## Catch the invisible result of ordihull...
pl <- ordihull(mod, Management, scaling = 3, label = TRUE)
## ... and find centres and areas of the hulls
summary(pl)
ordispider(mod, col="red", scaling = 3)
plot(mod, type = "p", display="sites")
ordicluster(mod, hclust(vegdist(dune)), prune=3, col = "blue")
plot(mod, type="n", display = "sites")
text(mod, display="sites", labels = as.character(Management))
pl <- ordiellipse(mod, Management, kind="se", conf=0.95, lwd=2, col="blue")
summary(pl)
}
\keyword{aplot }

