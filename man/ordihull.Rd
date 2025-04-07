\name{ordihull}
\alias{ordihull}
\alias{ordispider}
\alias{ordiellipse}
\alias{ordibar}
\alias{ordicluster}
\alias{summary.ordihull}
\alias{scores.ordihull}
\alias{summary.ordiellipse}
\alias{summary.ordiareatest}
\alias{ordiareatest}

\title{Display Groups or Factor Levels in Ordination Diagrams}

\description{ Functions to add convex hulls, \dQuote{spider} graphs, ellipses
  or cluster dendrogram to ordination diagrams. The ordination
  diagrams can be produced by \code{vegan} \code{\link{plot.cca}},
  \code{\link{plot.decorana}} or \code{\link{ordiplot}}.  }

\usage{
ordihull(ord, groups, display = "sites", draw = c("lines","polygon", "none"),
         col = NULL, alpha = 127, show.groups, label = FALSE,
         border = NULL, lty = NULL, lwd = NULL, ...)
ordiellipse(ord, groups, display="sites", kind = c("sd","se", "ehull"),
         conf, draw = c("lines","polygon", "none"),
	 w, col = NULL, alpha = 127, show.groups, label = FALSE,
         border = NULL, lty = NULL, lwd=NULL, ...)
ordibar(ord, groups, display = "sites", kind = c("sd", "se"), conf,
         w, col = 1, show.groups, label = FALSE, lwd = NULL, length = 0,  ...)
ordispider(ord, groups, display="sites", w, spiders = c("centroid", "median"),
         show.groups, label = FALSE, col = NULL, lty = NULL, lwd = NULL, ...)
ordicluster(ord, cluster, prune = 0, display = "sites",
            w, col = 1, draw = c("segments", "none"), ...)
\method{summary}{ordihull}(object, ...)
\method{summary}{ordiellipse}(object, ...)
ordiareatest(ord, groups, area = c("hull", "ellipse"), kind = "sd",
         permutations = 999, parallel = getOption("mc.cores"), ...)
\method{summary}{ordiareatest}(object, ...)
}

\arguments{
  \item{ord}{An ordination object or an \code{\link{ordiplot}} object. }
  \item{groups}{Factor giving the groups for which the graphical item is
    drawn. }
  \item{display}{Item to displayed. }

  \item{draw}{character; how should objects be represented on the plot?
    For \code{ordihull} and \code{ordiellipse} use either
    \code{\link{lines}} or \code{\link{polygon}} to draw the
    lines. For \code{ordicluster}, line segments are drawn using
    \code{\link{segments}}. To suppress plotting, use
    \code{"none"}. Graphical parameters are passed to both. The main
    difference is that \code{polygon}s may be filled and
    non-transparent. With \code{none} nothing is drawn, but the function
    returns the \code{\link{invisible}} plotting.}

  \item{col}{Colour of hull or ellipse lines (if \code{draw = "lines"})
    or their fills (if \code{draw = "polygon"}) in \code{ordihull} and
    \code{ordiellipse}.  When \code{draw = "polygon"}, the colour of
    bordering lines can be set with argument \code{border} of the
    \code{\link{polygon}} function. For other functions the effect
    depends on the underlining functions this argument is passed to.
    When multiple values of \code{col} are specified these are used
    for each element of \code{names(table(groups))} (in that order),
    shorter vectors are recycled. Function \code{ordicluster} has
    no \code{groups}, and there the argument will be recycled for
    points, and the colour of connecting lines is a mixture of point
    s in the cluster.}

  \item{alpha}{Transparency of the fill \code{col}our with \code{draw
    = "polygon"} in \code{ordihull} and \code{ordiellipse}.  The
    argument takes precedence over possible transparency definitions
    of the colour. The value must be in range \eqn{0...255}, and low
    values are more transparent.  Transparency is not available in all
    graphics devices or file formats.}

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

  \item{kind}{Draw standard deviations of points (\code{sd}), standard
    errors (\code{se}) or ellipsoid hulls that enclose all points in
    the group (\code{ehull}).}

  \item{conf}{Confidence limit for ellipses, e.g. 0.95. If given, the
    corresponding \code{sd} or \code{se} is multiplied with the
    corresponding value found from the Chi-squared distribution with
    2df. }

  \item{spiders}{Are centres or spider bodies calculated either as
    centroids (averages) or spatial medians.}

  \item{cluster}{Result of hierarchic cluster analysis, such as
    \code{\link{hclust}} or \code{\link[cluster]{agnes}}.}

  \item{prune}{Number of upper level hierarchies removed from the
    dendrogram. If \code{prune} \eqn{>0}, dendrogram will be
    disconnected.}

  \item{object}{A result object from \code{ordihull}, \code{ordiellipse}
    or \code{ordiareatest}. The result is \code{\link{invisible}}, but
    it can be saved, and used for summaries (areas etc. of hulls and
    ellipses). }

  \item{area}{Evaluate the area of convex hulls of \code{ordihull}, or of
    ellipses of \code{ordiellipse}.}

  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}

  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}

  \item{lty, lwd, border}{Vectors of these parameters can be supplied
    and will be applied (if appropriate) for each element of
    \code{names(table(groups))} (in that order). Shorter vectors will be
    recycled.}

  \item{length}{Width (in inches) of the small (\dQuote{caps}) at the
    ends of the bar segment (passed to \code{\link{arrows}}).}

  \item{\dots}{Parameters passed to other functions. }
}

\details{

  Function \code{ordihull} draws \code{\link{lines}} or
  \code{\link{polygon}}s for the convex
  hulls found by function \code{\link{chull}} encircling
  the items in the groups.

  Function \code{ordiellipse} draws \code{\link{lines}} or
  \code{\link{polygon}}s for ellipses by \code{groups}. The function
  can either draw standard deviation of points (\code{kind="sd"}) or
  standard error of the (weighted) centroids (\code{kind="se"}), and
  the (weighted) correlation defines the direction of the principal
  axis of the ellipse. When \code{kind = "se"} is used together with
  argument \code{conf}, the ellipses will show the confidence regions
  for the locations of group centroids. With \code{kind="ehull"} the
  function draws an ellipse that encloses all points of a group using
  \code{\link[cluster]{ellipsoidhull}} (\pkg{cluster} package).

  Function \code{ordibar} draws crossed \dQuote{error bars} using
  either either standard deviation of point scores or standard error
  of the (weighted) average of scores. These are the principal axes of
  the corresponding \code{ordiellipse}, and are found by principal
  component analysis of the (weighted) covariance matrix.

  Functions \code{ordihull} and \code{ordiellipse} return invisibly an
  object that has a \code{summary} method that returns the coordinates
  of centroids and areas of the hulls or ellipses. Function
  \code{ordiareatest} studies the one-sided hypothesis that these
  areas are smaller than with randomized \code{groups}. Argument
  \code{kind} can be used to select the kind of ellipse, and has no
  effect with convex hulls.

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

  Function \code{ordihull} and \code{ordiellipse} return invisibly an
  object that has a \code{summary} method that returns the coordinates
  of centroids and areas of the hulls or ellipses. Function
  \code{ordiareatest} studies the one-sided hypothesis that these
  areas are smaller than with randomized \code{groups}, and its
  \code{summary} is based on \code{\link{summary.permustats}} with a
  summary of permutations.

}

\author{ Jari Oksanen }

\seealso{The functions pass parameters to basic graphical functions,
  and you may wish to change the default values in
  \code{\link{lines}}, \code{\link{segments}} and
  \code{\link{polygon}}. You can pass parameters to
  \code{\link{scores}} as well. Underlying functions for
  \code{ordihull} is \code{\link{chull}}. The underlying function for
  ellipsoid hulls in \code{ordiellipse} is
  \code{\link[cluster]{ellipsoidhull}}. }

\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ Management, dune.env)
plot(mod, type="n", scaling = "symmetric")
## Catch the invisible result of ordihull...
pl <- with(dune.env, ordihull(mod, Management,
                              scaling = "symmetric", label = TRUE))
## ... and find centres and areas of the hulls
summary(pl)
## use more colours and add ellipsoid hulls
plot(mod, type = "n")
pl <- with(dune.env, ordihull(mod, Management,
                              scaling = "symmetric", col = 1:4,
                              draw="polygon", label =TRUE))
with(dune.env, ordiellipse(mod, Management, scaling = "symmetric",
                           kind = "ehull", col = 1:4, lwd=3))
## ordispider to connect WA and LC scores
plot(mod, dis=c("wa","lc"), type="p")
ordispider(mod)
## Other types of plots
plot(mod, type = "p", display="sites")
cl <- hclust(vegdist(dune))
ordicluster(mod, cl, prune=3, col = cutree(cl, 4))
## confidence ellipse: location of the class centroids
plot(mod, type="n", display = "sites")
with(dune.env, text(mod, display="sites", labels = as.character(Management),
                    col=as.numeric(Management)))
pl <- with(dune.env, ordiellipse(mod, Management, kind="se", conf=0.95, lwd=2,
                                 draw = "polygon", col=1:4, border=1:4,
                                 alpha=63))
summary(pl)
## add confidence bars
with(dune.env, ordibar(mod, Management, kind="se", conf=0.95, lwd=2, col=1:4,
                       label=TRUE))
}
\keyword{aplot }
