\name{scores.hclust}
\alias{scores.hclust}

\title{
  Coordinates of Leaves and Internal Nodes in a hclust Tree
}

\description{ The function finds the coordinates that will be used for
  internal nodes and leaves when an \code{\link{hclust}} object is
  plotted. These help in annotating the plotted dendrogram.  }

\usage{
scores.hclust(x, display = "internal", ...)
}

\arguments{
  \item{x}{
  An \code{\link{hclust}} result object.
}
  \item{display}{
  Return \code{"internal"} nodes or \code{"terminal"} nodes (also
  called \code{"leaves"}.
}
  \item{\dots}{
  Other arguments passed to the function (ignored).
}
}

\details{

  The function returns the coordinates of nodes in an
  \code{\link{hclust}} plot as two-column matrix. First column called
  \code{x} gives the horizontal coordinates which for \eqn{n} terminal
  nodes (leaves) is an integer sequence \eqn{1..n}. The second column
  called \code{height} gives the merge value. For terminal nodes
  (leaves) this the value at which the item is merged to the tree, and
  in plots the labels can still hang below this level, as defined by
  the argument \code{hang} in \code{\link{plot.hclust}}.

  The function only works with \code{\link{hclust}} objects; it does
  not work with \code{\link{dendrogram}}.

}

\value{
  A two-column matrix of coordinates.
}

\author{
  Jari Oksanen.
}

\note{
  This function may be removed as useless.  
}


\seealso{
  \code{\link{hclust}}, \code{\link{plot.hclust}}.
}
\examples{
## Show values that were used in reordering a tree
data(mite, mite.env)
hc <- hclust(vegdist(mite))
hc <- with(mite.env, reorder(hc, WatrCont))
with(mite.env, plot(hc, labels=round(WatrCont), cex=0.7))
ordilabel(scores(hc), label=round(hc$value), cex=0.7)
}
\keyword{ multivariate }

