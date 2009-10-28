\name{pcnm}
\alias{pcnm}
\title{ Function to compute classical PCNM (Principal Coordinates of
  Neighbourhood Matrix). }
\description{
  This function computed classical PCNM by the principal coordinate
  analysis of a truncated distance matrix. 
}
\usage{
pcnm(matdist, threshold, w)
}

\arguments{
  \item{dis}{ A distance matrix. }
  \item{threshold}{ A threshold value or truncation distance. If
    missing, minimum distance giving connected network will be
    used. This is found as the longest distance in the minimum spanning
    tree of \code{dis}. }
  \item{w}{Prior weights for rows.}
}

\details{
  Principal Coordinates of Neighbourhood Matrix (PCNM) map distances
  between rows onto rectangular matrix on rows using a truncation
  threshold for long distances (Borcard & Legendre 2002). If original
  distances were Euclidean distances in two dimensions (like normal
  spatial distances), they could be mapped onto two dimensions if there
  is no truncation of distances. Because of truncation, there will be a
  higher number of principal coordinates. The selection of truncation
  distance has a huge influence on the PCNM vectors. The default is to
  use the longest distance to keep data connected. The distances above
  truncation threshold are given an arbitrary value of 4 times
  threshold.  For regular data, the first PCNM vectorsshow a wide scale
  variation and later PCNM vectors show smaller scale variation (Borcard
  & Legendre 2002), but for irregular data the intepretation is not as
  clear.

  The PCNM functions are used to express distances in rectangular form
  that is similar to normal explanatory variables used in, e.g.,
  constrained ordination (\code{\link{rda}}, \code{\link{cca}} and
  \code{\link{capscale}}) or univariate regression (\code{\link{lm}})
  together with environmental variables (row weights should be supplied
  with \code{\link{cca}}; see Examples). This is regarded as a more
  powerful method than forcing rectangular environmental data into
  distances and using them in partial mantel analysis
  (\code{\link{mantel.partial}}) together with geographic distances
  (Legendre et al. 2008, but see Tuomisto & Ruokolainen 2008).
  
  The function is based on \code{pcnm} function in Dray's unreleased
  \pkg{spacemakeR} package. The differences are that the current
  function usesr \code{\link[vegan]{spantree}} as an internal support
  function. The current function also can use prior weights for rows by
  using weighted metric scaling of \code{\link[vegan]{wcmscale}} of
  \pkg{vegan}. The use of row weights allows finding orthonormal PCNMs
  also for correspondence analysis (e.g., \code{\link[vegan]{cca}} or
  \code{\link[ade4]{cca}}).
}

\value{
  A list of three elements:
  \item{values }{Eigenvalues obtained by the principal coordinates
    analysis.} 
  \item{vectors }{Eigenvectors obtained by the principal coordinates
    analysis. They have been normalized to unit norm.}
 \item{threshold}{Truncation distance.} 
}
\references{
  Borcard D. and Legendre P. (2002) All-scale spatial analysis of
  ecological data by means of principal coordinates of neighbour
  matrices. \emph{Ecological Modelling} \bold{153}, 51--68.
  
  Legendre, P., Bordard, D and Peres-Neto, P. (2008) Analyzing or
  explaining beta diversity? Comment. \emph{Ecology} \bold{89},
  3238--3244.

  Tuomisto, H. & Ruokolainen, K. (2008) Analyzing or explaining beta
  diversity? A reply. \emph{Ecology} \bold{89}, 3244--3256.
}

\author{Jari Oksanen, based on the code of Stephane Dray.}
\seealso{ \code{\link[vegan]{spantree}}. }
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
\keyword{ multivariate }
