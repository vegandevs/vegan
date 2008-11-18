\name{specpool}
\alias{specpool}
\alias{specpool2vect}
\alias{estimateR}
\alias{estimateR.default}
\alias{estimateR.matrix}
\alias{estimateR.data.frame}
\title{ Extrapolated Species Richness in a Species Pool}
\description{
  The functions estimate the extrapolated species richness in a species
  pool, or the number of unobserved species. Function \code{specpool}
  is based on incidences in sample sites, and gives a single estimate
  for a collection of sample sites (matrix).  Function \code{estimateR}
  is based on abundances (counts) on single sample site. 
}
\usage{
specpool(x, pool)
specpool2vect(X, index = c("Jack.1","Jack.2", "Chao", "Boot","Species"))
estimateR(x, ...)
}

\arguments{
  \item{x}{Data frame or matrix with species data.}
  \item{pool}{A vector giving a classification for pooling the sites in
    the species data. If missing, all sites are pooled together.}
  \item{X}{A \code{specpool} result object.}
  \item{index}{The selected index of extrapolated richness.}
  \item{...}{Other parameters (not used).}
}
\details{
  Many species will always remain unseen or undetected in a collection
  of sample plots.  The function uses some popular ways of estimating
  the number of these unseen species and adding them to the observed
  species richness (Palmer 1990, Colwell & Coddington 1994).

  The incidence-based estimates in \code{specpool} use the frequencies
  of species in a collection of sites.
  In the following, \eqn{S_P} is the extrapolated richness in a pool,
  \eqn{S_0} is the observed number of species in the
  collection, \eqn{a_1}{a1} and \eqn{a_2}{a2} are the number of species
  occurring only in one or only in two sites in the collection, \eqn{p_i}
  is the frequency of species \eqn{i}, and \eqn{N} is the number of
  sites in the collection.  The variants of extrapolated richness in
  \code{specpool} are:
  \tabular{ll}{
    Chao
    \tab \eqn{S_P = S_0 + a1^2/(2*a2)}
    \cr
    First order jackknife
    \tab \eqn{S_P = S_0 + a_1 \frac{N-1}{N}}{S_P = S_0 + a1*(N-1)/N}
    \cr
    Second order jackknife
    \tab \eqn{S_P = S_0 + a_1 \frac{2N - 3}{N} - a_2 \frac{(N-2)^2}{N
	(N-1)}}{S_P = S_0 + a1*(2*n-3)/n - a2*(n-2)^2/n/(n-1)}
    \cr
    Bootstrap
    \tab \eqn{S_P = S_0 + \sum_{i=1}^{S_0} (1 - p_i)^N}{S_P = S_0 + Sum
      (1-p_i)^N}
    }

    The abundance-based estimates in \code{estimateR} use counts (frequencies) of
    species in a single site. If called for a matrix or data frame, the
    function will give separate estimates for each site.  The two
    variants of extrapolated richness in \code{estimateR} are Chao
    (unbiased variant) and ACE.  In the Chao estimate
    \eqn{a_i} refers to number of species with abundance \eqn{i} instead
    of incidence: 
    \tabular{ll}{
    Chao
    \tab \eqn{S_P = S_0 + \frac{a_1 (a_1 -1)}{2 (a_2 + 1)}}{S_P = S_0 +
      a1*(a1-1)/(2*(a2+1))}
    \cr
    ACE
    \tab \eqn{S_P = S_{abund} + \frac{S_{rare}}{C_{ace}}+ \frac{a_1}{C_{ace}}
      \gamma^2_{ace}}{S_P = S_abund + S_rare/C_ace + a1/C_ace * gamma^2}
    \cr
    where \tab
    \eqn{C_{ace} = 1 - \frac{a_1}{N_{rare}}}{C_{ace} = 1- a1/N_{rare}}
    \cr
    \tab \eqn{\gamma^2_{ace} = \max \left[ \frac{S_{rare} \sum_{i=1}^{10}
      i(i-1)a_i}{C_{ace} N_{rare} (N_{rare} - 1)}-1, 0 \right]}{gamma^2 = \max(
    S_rare/C_ace (sum[i=1..10] i*(i-1)*a_i) / N_rare/(N_rare-1) -1 , 0)}
    }
    Here \eqn{a_i} refers to number of species with abundance \eqn{i}
    and  \eqn{S_{rare}}{S_rare} is the number of rare
    species, 
    \eqn{S_{abund}}{S_abund} is the number of abundant species, with an
    arbitrary 
    threshold of abundance 10 for rare species, and \eqn{N_{rare}}{N_rare} is
    the number 
    of individuals in rare species.

    Functions estimate the standard errors of the estimates. These
    only concern the number of added species, and assume that there is
    no variance in the observed richness.
    The equations of standard errors are too complicated to be reproduced in
    this help page, but they can be studied in the \R source code of the
    function.
    The standard error are based on the following sources: Chao (1987)
    for the Chao estimate and Smith and van Belle (1984) for the
    first-order Jackknife and the bootstrap (second-order jackknife is
    still missing). 
    The variance estimator of \eqn{S_{ace}}{S_ace} was
    developed by Bob O'Hara (unpublished).
    
}
\value{
  Function \code{specpool} returns a data frame with entries for
 observed richness
  and each of the indices for each class in \code{pool} vector.  The
  utility function \code{specpool2vect} maps the pooled values into
  a vector giving the value of selected \code{index} for each original
  site. Function \code{estimateR} returns the estimates and their
  standard errors for each site. 
}
\references{
  Chao, A. (1987). Estimating the population size for capture-recapture
  data with unequal catchability. \emph{Biometrics} 43, 783--791.
  
  Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial
  biodiversity through
  extrapolation. \emph{Phil. Trans. Roy. Soc. London} B 345, 101--118.

  Palmer, M.W. (1990). The estimation of species richness by
  extrapolation. \emph{Ecology} 71, 1195--1198.

  Smith, E.P & van Belle, G. (1984). Nonparametric estimation of
  species richness. \emph{Biometrics} 40, 119--129.
}
\author{Bob O'Hara (\code{estimateR}) and Jari Oksanen (\code{specpool}).}
\note{
  The functions are based on assumption that there is a species pool:
  The community is closed so that there is a fixed pool size \eqn{S_P}.
  Such cases may exist, although I have not seen them yet.  All indices
  are biased for open communities.

  See \url{http://viceroy.eeb.uconn.edu/EstimateS} for a more complete
  (and positive) discussion and alternative software for some platforms.
}
\seealso{\code{\link{veiledspec}}, \code{\link{diversity}}, \code{\link{beals}}. }
\examples{
data(dune)
data(dune.env)
attach(dune.env)
pool <- specpool(dune, Management)
pool
op <- par(mfrow=c(1,2))
boxplot(specnumber(dune) ~ Management, col="hotpink", border="cyan3",
 notch=TRUE)
boxplot(specnumber(dune)/specpool2vect(pool) ~ Management, col="hotpink",
 border="cyan3", notch=TRUE)
par(op)
data(BCI)
estimateR(BCI[1:5,])
}
\keyword{ univar }

