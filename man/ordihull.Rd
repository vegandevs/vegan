\name{ordihull}
\alias{ordihull}
\alias{ordispider}
\alias{ordiellipse}
\alias{ordicluster}
\alias{weights.cca}
\alias{weights.rda}
\alias{weights.decorana}
\alias{summary.ordihull}
\alias{scores.ordihull}
\alias{summary.ordiellipse}

\title{Display Groups or Factor Levels in Ordination Diagrams}

\description{ Functions to add convex hulls, `spider' graphs, ellipses
  or cluster dendrogram to ordination diagrams. The ordination
  diagrams can be produced by \code{vegan} \code{\link{plot.cca}},
  \code{\link{plot.decorana}} or \code{\link{ordiplot}}.  }

\usage{
ordihull(ord, groups, display = "sites", draw = c("lines","polygon", "none"),
         show.groups, label = FALSE,  ...)
ordiellipse(ord, groups, display="sites", kind = c("sd","se"), conf,
         draw = c("lines","polygon", "none"), w = weights(ord, display),
         col = NULL, show.groups, label = FALSE, ...)
ordispider(ord, groups, display="sites", w = weights(ord, display),
         show.groups, label = FALSE, ...)
ordicluster(ord, cluster, prune = 0, display = "sites",
         w = weights(ord, display), ...)
\method{summary}{ordihull}(object, ...)
\method{summary}{ordiellipse}(object, ...)
}

\arguments{
  \item{ord}{An ordination object or an \code{\link{ordiplot}} object. }
  \item{groups}{Factor giving the groups for which the graphical item is
    drawn. }
  \item{display}{Item to displayed. }

  \item{draw}{Use either \code{\link{lines}} or \code{\link{polygon}}
    to draw the lines. Graphical parameters are passed to both. The
    main difference is that \code{polygon}s may be filled and
    non-transparent. With \code{none} nothing is drawn, but the
    function returns the \code{\link{invisible}} plotting data.}

  \item{show.groups}{Show only given groups. This can be a vector, or
    \code{TRUE} if you want to show items for which condition is
    \code{TRUE}. This argument makes it possible to use different
    colours and line types for groups. The default is to show all groups. }

  \item{label}{Label the \code{groups} by their names in the centroid
    of the object. \code{ordiellipse} and \code{ordihull} use standard
    \code{\link{text}}, and others use \code{\link{ordilabel}}.}

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

  \item{col}{Colour of ellipses or ellipse fills in
    \code{ordiellipse}. For other functions the effect depends on the
    underlining functions this argument is passed to.}

  \item{prune}{Number of upper level hierarchies removed from the
    dendrogram. If \code{prune} \eqn{>0}, dendrogram will be
    disconnected.}

  \item{object}{A result object from \code{ordihull} or
    \code{ordiellipse}. The result is \code{\link{invisible}}, but it
    can be saved, and used for summaries (areas etc. of hulls and
    ellipses). }

  \item{\dots}{Parameters passed to graphical functions or to
    \code{\link{scores}} to select axes and scaling etc. } 
}

\details{
  Function \code{ordihull} draws \code{\link{lines}} or
  \code{\link{polygon}}s for the convex
  hulls found by function \code{\link{chull}} encircling
  the items in the groups. 
  
  Function \code{ordiellipse} draws \code{\link{lines}} or
  \code{\link{polygon}}s for dispersion ellipses
  using either standard deviation of point scores or standard error of
  the (weighted) average of scores, and the (weighted) correlation
  defines the direction of the principal axis of the ellipse. 
  An ellipsoid hull can be drawn with function
  \code{\link[cluster]{ellipsoidhull}} of package \pkg{cluster}.

  Function \code{ordispider} draws a \sQuote{spider} diagram where
  each point is connected to the group centroid with
  \code{\link{segments}}.  Weighted centroids are used in the
  correspondence analysis methods \code{\link{cca}} and
  \code{\link{decorana}} or if the user gives the weights in the
  call. If \code{ordispider} is called with \code{\link{cca}} or
  \code{\link{rda}} result without \code{groups} argument, the
  function connects each \sQuote{WA} scores to the corresponding
  \sQuote{LC} score. If the argument is a (\code{invisible})
  \code{ordihull} object, the function will connect the points of the
  hull to their centroid.

  Function \code{ordicluster} overlays a cluster dendrogram onto
  ordination. It needs the result from a hierarchic clustering such as
  \code{\link{hclust}} or \code{\link[cluster]{agnes}}, or other with
  a similar structure. Function \code{ordicluster} connects cluster
  centroids to each other with line \code{\link{segments}}. Function
  uses centroids of all points in the clusters, and is therefore
  similar to average linkage methods.  
}

\note{These functions add graphical items to ordination graph: You
  must draw a graph first. To draw line segments, grids or arrows, see
  \code{\link{ordisegments}}, \code{\link{ordigrid}}
  and\code{\link{ordiarrows}}.  }

\value{

  Functions \code{ordihull}, \code{ordiellipse} and \code{ordispider}
  return the \code{\link{invisible}} plotting structure. 

  Function \code{ordispider} return the coordinates to which each
  point is connected (centroids or \sQuote{LC} scores).

  Function \code{ordihull} returns a list of coordinates of the hulls
  (which can be extracted with \code{scores}), and \code{ordiellipse}
  returns a list of covariance matrices and scales used in drawing the
  ellipses.  These result objects have a \code{summary} method that
  returns the coordinates of the centres of the ellipses or hulls and
  their surface areas in user units. The centres of the hulls may
  differ from the location of the \code{label} which is the centre of
  the points instead of the centre of the polygon. With \code{draw =
  "none"} only the result object is returned and nothing is drawn.

}

\author{ Jari Oksanen }

\seealso{The functions pass parameters to basic graphical functions,
  and you may wish to change the default values in
  \code{\link{lines}}, \code{\link{segments}} and
  \code{\link{polygon}}. You can pass parameters to
  \code{\link{scores}} as well. Underlying function for
  \code{ordihull} is \code{\link{chull}}.  }

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
## use ordispider to label and mark the hull
plot(mod, type = "n")
pl <- ordihull(mod, Management, scaling = 3)
ordispider(pl, col="red", lty=3, label = TRUE )
## ordispider to connect WA and LC scores
plot(mod, dis=c("wa","lc"), type="p")
ordispider(mod)
## Other types of plots
plot(mod, type = "p", display="sites")
ordicluster(mod, hclust(vegdist(dune)), prune=3, col = "blue")
plot(mod, type="n", display = "sites")
text(mod, display="sites", labels = as.character(Management))
pl <- ordiellipse(mod, Management, kind="se", conf=0.95, lwd=2, col="blue")
summary(pl)
}
\keyword{aplot }
