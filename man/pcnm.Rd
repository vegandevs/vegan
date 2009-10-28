\name{pcnm}
\alias{pcnm}
\title{ Function to compute classical PCNM(Principal Coordinates of Neighborhood Matrix). }
\description{
  This function computed classical PCNM by the principal coordinate analysis of a truncated distance matrix.
}
\usage{
pcnm(matdist, threshold, support = c("vegan", "ade4"))
}

\arguments{
  \item{matdist}{ A distance matrix. }
  \item{threshold}{ A threshold value or truncation distance. If missing,
    minimum distance giving connected network will be used. This is found
    as the longest distance in the minimum spanning tre of
    \code{matdist}. }
  \item{support}{Use either \pkg{vegan} or \pkg{ade4} function for the
    minimum spanning tree.}
}

\value{
  A list of three elements:
  \item{values }{Eigenvalues obtained by the principal coordinates analysis.}
  \item{vectors }{Eigenvectors obtained by the principal coordinates analysis. They have been normalized to unit norm.}
 \item{threshold}{Truncation distance.} 
}
\references{ Borcard D. and Legendre P. (2002) All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. \emph{Ecological Modelling} \bold{153}, 51--68. }
\author{ Stephane Dray, adapted to vegan by Jari Oksanen }
\seealso{ \code{\link[vegan]{spantree}}, \code{\link[ade4]{mstree}}. }
\examples{
data(mite.xy)
pcnm1 <- pcnm(dist(mite.xy))
op <- par(mfrow=c(1,3))
ordisurf(mite.xy, pcnm1$vectors[,1])
ordisurf(mite.xy, pcnm1$vectors[,2])
ordisurf(mite.xy, pcnm1$vectors[,3])
par(op)

}
\keyword{ spatial }
