\name{spantree}
\alias{spantree}
\alias{cophenetic.spantree}
\alias{plot.spantree}
\alias{lines.spantree}

\title{Minimum Spanning Tree}
\description{
  Function \code{spantree} finds a minimum spanning tree
  connecting all points, but disregarding dissimilarities that are at or
  above the threshold or \code{NA}.  
}
\usage{
spantree(dis, toolong = 0)
\method{cophenetic}{spantree}(x)
\method{plot}{spantree}(x, ord, cex = 0.7, type = "p", labels, dlim,
     FUN = sammon,  ...)
\method{lines}{spantree}(x, ord, display="sites", ...)
}

\arguments{
  \item{dis}{Dissimilarity data inheriting from class \code{dist} or
    a an object, such as a matrix, that can be converted to a
    dissimilarity matrix. Functions \code{\link{vegdist}} and
    \code{\link{dist}} are some functions producing suitable
    dissimilarity data.}
  \item{toolong}{ Shortest dissimilarity regarded as \code{NA}.
    The function uses a fuzz factor, so
    that dissimilarities close to the limit will be made \code{NA}, too.
    If \code{toolong = 0} (or negative), no dissimmilarity is regarded
    as too long.
  }
  \item{x}{A \code{spantree} result object.}
  \item{ord}{An ordination configuration, or an ordination result known
    by \code{\link{scores}}.}
  \item{cex}{Character expansion factor.}
  \item{type}{Observations are plotted as points with
    \code{type="p"} or \code{type="b"}, or as text label with
    \code{type="t"}. The tree (lines) will always be plotted.}
  \item{labels}{Text used with \code{type="t"} or node names if this is
    missing.}
  \item{dlim}{A ceiling value used to highest \code{cophenetic} dissimilarity.}
  \item{FUN}{Ordination function to find the configuration from
    cophenetic dissimilarities. }
  \item{display}{Type of \code{\link{scores}} used for \code{ord}.}
  \item{\dots}{Other parameters passed to functions.}
}
\details{
  
  Function \code{spantree} finds a minimum spanning tree for
  dissimilarities (there may be several minimum spanning trees, but the
  function finds only one). Dissimilarities at or above the threshold
  \code{toolong} and \code{NA}s are disregarded, and the spanning tree
  is found through other dissimilarities. If the data are disconnected,
  the function will return a disconnected tree (or a forest), and the
  corresponding link is \code{NA}. Connected subtrees can be identified
  using \code{\link{distconnected}}.

  Function \code{cophenetic} finds distances between all points along
  the tree segments. Function \code{plot} displays the tree over a
  supplied ordination configuration, and \code{lines} adds a spanning
  tree to an ordination graph. If configuration is not supplied for \code{plot},
  the function ordinates the cophenetic dissimilarities of the
  spanning tree and overlays the tree on this result. The default
  ordination function is \code{\link[MASS]{sammon}} (package \pkg{MASS}),
  because Sammon scaling emphasizes structure in the neighbourhood of
  nodes and may be able to beautifully represent the tree (you may need
  to set \code{dlim}, and sometimes the results will remain
  twisted). These ordination methods do not work with disconnected
  trees, but you must supply the ordination configuration. Function
  \code{lines} will overlay the tree in an existing plot. 

  Function \code{spantree} uses Prim's method
  implemented as priority-first search for dense graphs (Sedgewick
  1990). Function \code{cophenetic} uses function
  \code{\link{stepacross}} with option \code{path = "extended"}. The
  \code{spantree} is very fast, but \code{cophenetic} is slow in very
  large data sets. 
}
\value{
  Function \code{spantree}
  returns an object of class \code{spantree} which is a
  list with two vectors, each of length \eqn{n-1}. The
  number of links in a tree is one less the number of observations, and
  the first item is omitted. The items are  
  \item{kid }{The child node of the parent, starting from parent number
    two. If there is no link from the parent, value will be \code{NA}
    and tree is disconnected at the node.}
  \item{dist }{Corresponding distance. If \code{kid = NA}, then
    \code{dist = 0}.}
  \item{labels }{Names of nodes as found from the input dissimilarities.} 
}
\references{
 Sedgewick, R. (1990). \emph{Algorithms in C}. Addison Wesley. 
}
\author{ Jari Oksanen }
\note{
  In principle, minimum spanning tree is equivalent to single linkage
  clustering that can be performed using \code{\link{hclust}} or
  \code{\link[cluster]{agnes}}. However, these functions combine
  clusters to each other and the information of the actually connected points
  (the ``single link'') cannot be recovered from the result. The
  graphical output of a single linkage clustering plotted with
  \code{\link{ordicluster}} will look very different from an equivalent
  spanning tree plotted with \code{lines.spantree}.
}


\seealso{\code{\link{vegdist}} or \code{\link{dist}} for getting
  dissimilarities,  and \code{\link{hclust}} or
    \code{\link[cluster]{agnes}} for single linkage clustering. 
}
\examples{
data(dune)
dis <- vegdist(dune)
tr <- spantree(dis)
## Add tree to a metric scaling 
plot(tr, cmdscale(dis), type = "t")
## Find a configuration to display the tree neatly
plot(tr, type = "t")
}
\keyword{ multivariate}

