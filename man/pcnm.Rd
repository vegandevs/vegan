\name{pcnm}
\alias{pcnm}
\title{ Function to compute classical PCNM(Principal Coordinates of Neighborhood Matrix). }
\description{
  This function computed classical PCNM by the principal coordinate
  analysis of a truncated distance matrix. 
}
\usage{
pcnm(matdist, threshold, support = c("vegan", "ade4"), w)
}

\arguments{
  \item{matdist}{ A distance matrix. }
  \item{threshold}{ A threshold value or truncation distance. If missing,
    minimum distance giving connected network will be used. This is found
    as the longest distance in the minimum spanning tre of
    \code{matdist}. }
  \item{support}{Use either \pkg{vegan} or \pkg{ade4} function for the
    minimum spanning tree.}
  \item{w}{Prior weights for rows.}
}

\details{
  The function is based on \code{pcnm} function in Dray's
  \pkg{spacemakeR} package. The differences are that the current
  function can use either \code{\link[vegan]{spantree}} of \pkg{vegan}
  or \code{\link[ade4]{mstree}} and \code{\link[ade4]{neig2mat}} of
  \pkg{ade4} as internal support function. The current function also can
  use prior weights for rows by using weighted metric scaling of
  \code{\link[vegan]{wcmscale}} of \pkg{vegan}. The use of row weights
  allows finding orthonormal PCNMs also for correspondence analysis
  (e.g., \code{\link[vegan]{cca}} or \code{\link[ade4]{cca}}).  }

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
## Weighted PCNM for CCA
data(mite)
rs <- rowSums(mite)/sum(mite)
pcnmw <- pcnm(dist(mite.xy), w = rs)
ord <- cca(mite ~ pcnmw$vectors)
}
\keyword{ spatial }
