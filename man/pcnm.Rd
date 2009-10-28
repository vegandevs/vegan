\name{pcnm}
\alias{pcnm}
\title{ Function to compute classical PCNM(Principal Coordinates of Neighborhood Matrix). }
\description{
  This function computed classical PCNM by the principal coordinate analysis of a truncated distance matrix.
}
\usage{
pcnm(matdist, thresh = give.thresh(as.dist(matdist)))
}

\arguments{
  \item{matdist}{ A distance matrix. }
  \item{thresh}{ A threshold value. }
}

\value{
  A list of two elements:
  \item{values }{Eigenvalues obtained by the principal coordinates analysis.}
  \item{vectors }{Eigenvectors obtained by the principal coordinates analysis. They have been normalized to unit norm.}
  
}
\references{ Borcard D. and Legendre P. (2002) All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. \emph{Ecological Modelling} \bold{153}, 51--68. }
\author{ Stephane Dray }
\seealso{ \code{\link{give.thresh}} }
\examples{
data(oribatid)
pcnm1 <- pcnm(dist(oribatid$xy))
par(mfrow=c(1,3))
s.value(oribatid$xy,pcnm1$vectors[,1],cleg=0, sub="PCNM 1",csub=3)
s.value(oribatid$xy,pcnm1$vectors[,2],cleg=0, sub="PCNM 2",csub=3)
s.value(oribatid$xy,pcnm1$vectors[,3],cleg=0, sub="PCNM 3",csub=3)


}
\keyword{ spatial }
