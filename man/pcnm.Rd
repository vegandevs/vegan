\name{pcnm}
\alias{pcnm}
\alias{scores.pcnm}
\title{ Principal Coordinates of Neighbourhood Matrix }
\description{
  This function computed classical PCNM by the principal coordinate
  analysis of a truncated distance matrix. These are commonly used to
  transform (spatial) distances to rectangular data that suitable for
  constrained ordination or regression. 
}
\usage{
pcnm(dis, threshold, w, dist.ret = FALSE)
}

\arguments{
  \item{dis}{ A distance matrix. }
  \item{threshold}{ A threshold value or truncation distance. If
    missing, minimum distance giving connected network will be
    used. This is found as the longest distance in the minimum spanning
    tree of \code{dis}. }
  \item{w}{Prior weights for rows.}
  \item{dist.ret}{Return the distances used to calculate the PCNMs.}
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
  function usesr \code{\link{spantree}} as an internal support
  function. The current function also can use prior weights for rows by
  using weighted metric scaling of \code{\link{wcmdscale}}. The use of
  row weights allows finding orthonormal PCNMs also for correspondence
  analysis (e.g., \code{\link{cca}}).
  }

\value{
  A list of the following elements:
  \item{values }{Eigenvalues obtained by the principal coordinates
    analysis.} 
  \item{vectors }{Eigenvectors obtained by the principal coordinates
    analysis. They are scaled to unit norm. The vectors can be extracted 
    with \code{scores} function. The default is to return all PCNM vectors,
    but argument \code{choices} selects the given vectors.} 
 \item{threshold}{Truncation distance.}
 \item{dist}{The distance matrix where values above \code{threshold}
    are replaced with arbitrary value of four times the
    threshold. String \code{"pcnm"} is added to the \code{method}
    attribute, and new attribute \code{threshold} is added to the
    distances. This is returned only when \code{dist.ret = TRUE}.  }
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
## Example from Borcard & Legendre (2002)
data(mite.xy)
pcnm1 <- pcnm(dist(mite.xy))
op <- par(mfrow=c(1,3))
## Map of PCNMs in the sample plot
ordisurf(mite.xy, scores(pcnm1, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(mite.xy, scores(pcnm1, choi=2), bubble = 4, main = "PCNM 2")
ordisurf(mite.xy, scores(pcnm1, choi=3), bubble = 4, main = "PCNM 3")
par(op)
## Plot first PCNMs against each other
ordisplom(pcnm1, choices=1:4)
## Weighted PCNM for CCA
data(mite)
rs <- rowSums(mite)/sum(mite)
pcnmw <- pcnm(dist(mite.xy), w = rs)
ord <- cca(mite ~ scores(pcnmw))
## Multiscale ordination: residual variance should have no distance
## trend
msoplot(mso(ord, mite.xy))
}
\keyword{ spatial }
\keyword{ multivariate }
