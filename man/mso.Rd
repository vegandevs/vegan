\name{mso}
\alias{mso}
\alias{print.mso}
\alias{msoplot}

\title{ Functions for performing and displaying a spatial partitioning
of cca or rda results}

\description{ The function \code{mso} adds an attribute \code{vario} to
  an object of class \code{"cca"} that describes the spatial
  partitioning of the \code{\link{cca}} object and performs an optional
  permutation test for the spatial independence of residuals. The
  function \code{plot.mso} creates a diagnostic plot of the spatial
  partitioning of the \code{"cca"} object.  }

\usage{
mso(object.cca, object.xy, grain = 1, round.up = FALSE, permutations = FALSE)
msoplot(x, alpha = 0.05, explained = FALSE, ...)
}
\arguments{
  \item{object.cca}{ An object of class cca, created by the \code{\link{cca}} or
    \code{\link{rda}} function.} 
  \item{object.xy}{ A vector, matrix or data frame with the spatial
    coordinates of the data represented by object.cca. Must have the
    same number of rows as \code{object.cca\$CA\$Xbar} (see
    \code{\link{cca.object}}). } 
  \item{grain}{ Interval size for distance classes.}
  \item{round.up}{ Determines the choice of breaks. If false, distances
    are rounded to the nearest multiple of grain. If true, distances are
    rounded to the upper multiple of grain.} 
  \item{permutations}{ If false, suppresses the permutation test. If an
    integer, determines the number of permutations for the Mantel test
    of spatial independence of residual inertia.}
  \item{x}{A result object of \code{mso}.}
  \item{alpha}{ Significance level for the two-sided permutation test of
    the Mantel statistic for spatial independence of residual inertia
    and for the point-wise envelope of the variogram of the total
    variance. A Bonferroni-type correction can be achieved by dividing
    the overall significance value (e.g. 0.05) by the number of distance
    classes.} 
  \item{explained}{ If false, suppresses the plotting of the variogram
    of explained variance.}
  \item{\dots}{Other arguments passed to functions.}
}
\details{
  The Mantel test is an adaptation of the function \code{\link{mantel}} of the
  \pkg{vegan} package to the parallel testing of several distance classes. It
  compares the mean inertia in each distance class to the pooled mean
  inertia of all other distance classes. 

  If there are explanatory variables (RDA, CCA, pRDA, pCCA) and a
  significance test for residual autocorrelation was performed when
  running the function \code{mso}, the function \code{plot.mso} will
  print an estimate of how much the autocorrelation (based on
  significant distance classes) causes the global error variance of the
  regression analysis to be underestimated 

}
\value{
  The function \code{mso} returns an amended \code{cca} or \code{rda}
  object with the additional attributes \code{grain}, \code{H},
  \code{H.test} and \code{vario}.  
  \item{grain }{ The grain attribute defines the interval size of the
    distance classes .} 
  \item{H }{ H is an object of class 'dist' and contains the geographic
    distances between observations.} 
  \item{H.test }{ H.test contains a set of dummy variables that describe
    which pairs of observations (rows = elements of \code{object\$H}) fall in
    which distance class (columns). } 
  \item{vario }{ The vario attribute is a data frame that contains some
    or all of the following components for the rda case (cca case in
    brackets): 
    \item{H }{ Distance class as multiples of grain.}
    \item{Dist }{ Average distance of pairs of observations in distance class H.}
    \item{n }{ Number of unique pairs of observations in distance class
      H.} 
    \item{All }{ Empirical (chi-square) variogram of total variance
      (inertia).} 
    \item{Sum }{ Sum of empirical (chi-square) variograms of explained
      and residual variance (inertia).} 
    \item{CA }{ Empirical (chi-square) variogram of residual variance
      (inertia).} 
    \item{CCA }{ Empirical (chi-square) variogram of explained variance
      (inertia).} 
    \item{pCCA }{ Empirical (chi-square) variogram of conditioned
      variance (inertia).} 
    \item{se }{ Standard error of the empirical (chi-square) variogram
      of total variance (inertia).} 
    \item{CA.signif }{P-value of permutation test for spatial
      independence of residual variance (inertia).} 
  }
}
\references{ Wagner, H.H. 2004. Direct multi-scale ordination with
  canonical correspondence analysis. \emph{Ecology} 85: 342--351. } 
\author{ The responsible author was Helene Wagner.} 
\note{ The function is based on the code published in the Ecological
  Archives E085-006 (\url{http://www.esapubs.org/archive/ecol/E085/006/default.htm}).  }
\seealso{ Function \code{\link{cca}} and \code{\link{rda}},
  \code{\link{cca.object}}. }
\examples{
## Reconstruct worked example of Wagner (submitted):
X <- matrix(c(1, 2, 3, 2, 1, 0), 3, 2)
Y <- c(3, -1, -2)
tmat <- c(1:3)
## Canonical correspondence analysis (cca):
Example.cca <- cca(X, Y)
Example.cca <- mso(Example.cca, tmat)
msoplot(Example.cca)
Example.cca$vario

## Correspondence analysis (ca):
Example.ca <- mso(cca(X), tmat)
msoplot(Example.ca)

## Unconstrained ordination with test for autocorrelation
## using oribatid mite data set as in Wagner (2004)
data(mite)
data(mite.env)
data(mite.xy)

mite.cca <- cca(log(mite + 1))
mite.cca <- mso(mite.cca, mite.xy, grain =  1, permutations = 100)
msoplot(mite.cca)
mite.cca

## Constrained ordination with test for residual autocorrelation
## and scale-invariance of species-environment relationships
mite.cca <- cca(log(mite + 1) ~ SubsDens + WatrCont + Substrate + Shrub + Topo, mite.env)
mite.cca <- mso(mite.cca, mite.xy, permutations = 100)
msoplot(mite.cca)
mite.cca
}
\keyword{ spatial }
\keyword{ multivariate }
