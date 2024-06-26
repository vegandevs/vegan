\name{mso}
\alias{mso}
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
mso(object.cca, object.xy, grain = 1, round.up = FALSE, permutations = 0)
msoplot(x, alpha = 0.05, explained = FALSE, ylim = NULL, legend = "topleft", ...)
}
\arguments{
  \item{object.cca}{ An object of class cca, created by the \code{\link{cca}} or
    \code{\link{rda}} function.}
   \item{object.xy}{ A vector, matrix or data frame with the spatial
    coordinates of the data represented by \code{object.cca}. The
    number of rows must match the number of observations (as given by
    \code{nobs}) in \code{cca.object}. Alternatively, interpoint
    distances can be supplied as a \code{\link{dist}} object. }
   \item{grain}{ Interval size for distance classes.}
  \item{round.up}{ Determines the choice of breaks. If false, distances
    are rounded to the nearest multiple of grain. If true, distances are
    rounded to the upper multiple of grain.} 
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{x}{A result object of \code{mso}.}
  \item{alpha}{ Significance level for the two-sided permutation test of
    the Mantel statistic for spatial independence of residual inertia
    and for the point-wise envelope of the variogram of the total
    variance. A Bonferroni-type correction can be achieved by dividing
    the overall significance value (e.g. 0.05) by the number of distance
    classes.} 
  \item{explained}{ If false, suppresses the plotting of the variogram
    of explained variance.}
  \item{ylim}{Limits for y-axis.}
  \item{legend}{The x and y co-ordinates to be used to position the legend. 
    They can be specified by keyword or in any way which is accepted 
    by \code{\link{legend}}.}
  \item{\dots}{Other arguments passed to functions.}
}

\details{

  The Mantel test is an adaptation of the function
  \code{\link{mantel}} to the parallel testing of several distance
  classes and similar to multivariate \code{\link{mantel.correlog}}.
  It compares the mean inertia in each distance class to the pooled
  mean inertia of all other distance classes.

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
  \item{grain}{The grain attribute defines the interval size of the
    distance classes .} 
  \item{H }{ H is an object of class 'dist' and contains the geographic
    distances between observations.} 
  \item{H.test }{ H.test contains a set of dummy variables that describe
    which pairs of observations (rows = elements of \code{object$H}) fall in
    which distance class (columns). } 
  \item{vario }{ The vario attribute is a data frame that contains some
    or all of the following components for the rda case (cca case in
      brackets):
    \describe{
      \item{\code{H}}{Distance class as multiples of grain.}
      \item{\code{Dist}}{ Average distance of pairs of observations in distance class H.}
      \item{n }{ Number of unique pairs of observations in distance class
	H.} 
      \item{\code{All}}{ Empirical (chi-square) variogram of total variance
	(inertia).} 
      \item{\code{Sum}}{ Sum of empirical (chi-square) variograms of explained
	and residual variance (inertia).} 
      \item{\code{CA}}{ Empirical (chi-square) variogram of residual variance
	(inertia).} 
      \item{\code{CCA}}{ Empirical (chi-square) variogram of explained variance
	(inertia).} 
      \item{\code{pCCA}}{ Empirical (chi-square) variogram of conditioned
	variance (inertia).} 
      \item{\code{se}}{ Standard error of the empirical (chi-square) variogram
	of total variance (inertia).} 
      \item{\code{CA.signif}}{P-value of permutation test for spatial
	independence of residual variance (inertia).} 
    }
  }
}
  
\references{ Wagner, H.H. 2004. Direct multi-scale ordination with
  canonical correspondence analysis. \emph{Ecology} 85: 342--351. } 
\author{ The responsible author was Helene Wagner.} 
\note{ The function is based on the code published in the Ecological
  Archives E085-006 (\doi{10.1890/02-0738}).  }

\seealso{ \code{\link{mantel.correlog}}. }

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
mite.cca <- mso(mite.cca, mite.xy, grain =  1, permutations = 99)
msoplot(mite.cca)
mite.cca

## Constrained ordination with test for residual autocorrelation
## and scale-invariance of species-environment relationships
mite.cca <- cca(log(mite + 1) ~ SubsDens + WatrCont + Substrate + Shrub + Topo, mite.env)
mite.cca <- mso(mite.cca, mite.xy, permutations = 99)
msoplot(mite.cca)
mite.cca
}
\keyword{ spatial }
\keyword{ multivariate }
