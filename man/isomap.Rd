\name{isomap}
\alias{isomap}
\alias{isomapdist}
\alias{plot.isomap}
\alias{summary.isomap}

\title{ Isometric Feature Mapping Ordination }
\description{
  The function performs isometric feature mapping which consists of
  three simple steps: (1) retain only some of the shortest
  dissimilarities among objects, (2) estimate all dissimilarities as
  shortest path distances, and (3) perform metric scaling (Tenenbaum et
  al. 2000).
}
\usage{
isomap(dist, ndim=10, ...)
isomapdist(dist, epsilon, k, path = "shortest", fragmentedOK =FALSE, ...)
\method{summary}{isomap}(object, axes = 4, ...)
\method{plot}{isomap}(x, net = TRUE, n.col = "gray", type = "points", ...)
}

\arguments{
  \item{dist}{Dissimilarities. }
  \item{ndim}{Number of axes in metric scaling (argument \code{k} in 
   \code{\link{cmdscale}}).}
  \item{epsilon}{Shortest dissimilarity retained. }
  \item{k}{Number of shortest dissimilarities retained for a point. If
   both \code{epsilon} and \code{k} are given, \code{epsilon} will be used.  }
  \item{path}{Method used in \code{\link{stepacross}} to estimate the
    shortest path, with alternatives \code{"shortest"} and \code{"extended"}. }
  \item{fragmentedOK}{What to do if dissimilarity matrix is
    fragmented. If \code{TRUE}, analyse the largest connected group,
    otherwise stop with error. }
  \item{x, object}{An \code{isomap} result object.}
  \item{axes}{Number of axes displayed.}
  \item{net}{Draw the net of retained dissimilarities.}
  \item{n.col}{Colour of drawn net segments.}

  \item{type}{Plot observations either as \code{"points"},
    \code{"text"} or use \code{"none"} to plot no observations. The
    \code{"text"} will use \code{\link{ordilabel}} if \code{net = TRUE} 
    and \code{\link{ordiplot}} if \code{net = FALSE}, and pass
    extra arguments to these functions.}

  \item{\dots}{Other parameters passed to functions. }
}
\details{
  The function \code{isomap} first calls function \code{isomapdist} for
  dissimilarity transformation, and then performs metric scaling for the
  result. All arguments to \code{isomap} are passed to
  \code{isomapdist}. The functions are separate so that the
  \code{isompadist} transformation could be easily used with other
  functions than simple linear mapping of \code{\link{cmdscale}}.
  
  Function \code{isomapdist} retains either dissimilarities equal or shorter to
  \code{epsilon}, or if \code{epsilon} is not given, at least \code{k}
  shortest dissimilarities for a point.  Then a complete dissimilarity
  matrix is reconstructed using \code{\link{stepacross}} using either
  flexible shortest paths or extended dissimilarities (for details, see
  \code{\link{stepacross}}).

  De'ath (1999) actually published essentially the same method before
  Tenenbaum et al. (2000), and De'ath's function is available in  
  \code{\link[mvpart]{xdiss}} in package \pkg{mvpart}. The differences are that
  \code{isomap} introduced the \code{k} criterion, whereas De'ath only
  used \code{epsilon} criterion.  In practice, De'ath also retains
  higher proportion of dissimilarities than typical \code{isomap}.

  The \code{plot} function uses internally \code{\link{ordiplot}},
  except that it adds text over net using \code{\link{ordilabel}}. The
  \code{plot} function passes extra arguments to these functions.  In
  addition, \pkg{vegan3d} package has function
  \code{\link[vegan3d]{rgl.isomap}} to make dynamic 3D plots that can
  be rotated on the screen.
}

\value{
  Function \code{isomapdist} returns a dissimilarity object similar to
  \code{dist}. Function \code{isomap} returns an object of class
  \code{isomap} with \code{plot} and \code{summary} methods. The
  \code{plot} function returns invisibly an object of class
  \code{\link{ordiplot}}. Function \code{\link{scores}} can extract
  the ordination scores.
}
\references{
    De'ath, G. (1999)  Extended dissimilarity: a method of robust
     estimation of ecological distances from high beta diversity data.
     \emph{Plant Ecology} 144, 191--199

     Tenenbaum, J.B., de Silva, V. & Langford, J.C. (2000) A global
     network framework for nonlinear dimensionality
     reduction. \emph{Science} 290, 2319--2323.
}
\author{ Jari Oksanen }
\note{ 
  Tenenbaum et al. (2000) justify \code{isomap} as a tool of unfolding a
  manifold (e.g. a 'Swiss Roll'). Even with a manifold structure, the
  sampling must be even and dense so
  that dissimilarities along a manifold are shorter than across the
  folds. If data do not have such a manifold structure, the results are
  very sensitive to parameter values. 
}
\seealso{The underlying functions that do the proper work are
  \code{\link{stepacross}}, \code{\link{distconnected}} and \code{\link{cmdscale}}.  
  Package \pkg{mvpart} provides a parallel (but a bit different) implementation
  (\code{\link[mvpart]{xdiss}}).  Moreover, \pkg{vegan} function
  \code{\link{metaMDS}} may trigger \code{\link{stepacross}}
  transformation, but usually only  for longest dissimilarities.  The
  \code{plot} method of \pkg{vegan} minimum spanning tree function
  (\code{\link{spantree}}) has even more extreme way of isomapping things. }
\examples{
## The following examples also overlay minimum spanning tree to
## the graphics in red.
op <- par(mar=c(4,4,1,1)+0.2, mfrow=c(2,2))
data(BCI)
dis <- vegdist(BCI)
tr <- spantree(dis)
pl <- ordiplot(cmdscale(dis), main="cmdscale")
lines(tr, pl, col="red")
ord <- isomap(dis, k=3)
ord
pl <- plot(ord, main="isomap k=3")
lines(tr, pl, col="red")
pl <- plot(isomap(dis, k=5), main="isomap k=5")
lines(tr, pl, col="red")
pl <- plot(isomap(dis, epsilon=0.45), main="isomap epsilon=0.45")
lines(tr, pl, col="red")
par(op)
}
\keyword{ multivariate}

