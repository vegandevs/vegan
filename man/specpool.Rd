\name{specpool}
\alias{specpool}
\alias{specpool2vect}
\alias{poolaccum}
\alias{summary.poolaccum}
\alias{plot.poolaccum}
\alias{estimateR}
\alias{estimateR.default}
\alias{estimateR.matrix}
\alias{estimateR.data.frame}
\alias{estaccumR}

\title{ Extrapolated Species Richness in a Species Pool}
\description{
  The functions estimate the extrapolated species richness in a species
  pool, or the number of unobserved species. Function \code{specpool}
  is based on incidences in sample sites, and gives a single estimate
  for a collection of sample sites (matrix).  Function \code{estimateR}
  is based on abundances (counts) on single sample site. 
}
\usage{
specpool(x, pool, smallsample = TRUE)
estimateR(x, ...)
specpool2vect(X, index = c("jack1","jack2", "chao", "boot","Species"))
poolaccum(x, permutations = 100, minsize = 3)
estaccumR(x, permutations = 100, parallel = getOption("mc.cores"))
\method{summary}{poolaccum}(object, display, alpha = 0.05, ...)
\method{plot}{poolaccum}(x, alpha = 0.05, type = c("l","g"), ...)
}

\arguments{
  \item{x}{Data frame or matrix with species data or the analysis result 
    for \code{plot} function.}
  \item{pool}{A vector giving a classification for pooling the sites in
    the species data. If missing, all sites are pooled together.}
  \item{smallsample}{Use small sample correction \eqn{(N-1)/N}, where
    \eqn{N} is the number of sites within the \code{pool}.}
  \item{X, object}{A \code{specpool} result object.}
  \item{index}{The selected index of extrapolated richness.}
  \item{permutations}{Usually an integer giving the number
    permutations, but can also be a list of control values for the
    permutations as returned by the function \code{\link[permute]{how}}, 
    or a permutation matrix where each row gives the permuted indices.}
  \item{minsize}{Smallest number of sampling units reported.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}
  \item{display}{Indices to be displayed.}
  \item{alpha}{Level of quantiles shown. This proportion will be left outside
    symmetric limits.}
  \item{type}{Type of graph produced in \code{\link[lattice]{xyplot}}.}
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
    \tab \eqn{S_P = S_0 + \frac{a_1^2}{2 a_2}\frac{N-1}{N}}{S_P = S_0 + a1^2/(2*a2) * (N-1)/N}
    \cr
    Chao bias-corrected
    \tab \eqn{S_P = S_0 + \frac{a_1(a_1-1)}{2(a_2+1)} \frac{N-1}{N}}{S_P = S_0 + a1*(a1-1)/(2*(a2+1)) * (N-1)/N}
    \cr
    First order jackknife
    \tab \eqn{S_P = S_0 + a_1 \frac{N-1}{N}}{S_P = S_0 + a1*(N-1)/N}
    \cr
    Second order jackknife
    \tab \eqn{S_P = S_0 + a_1 \frac{2N - 3}{N} - a_2 \frac{(N-2)^2}{N
	(N-1)}}{S_P = S_0 + a1*(2*N-3)/N - a2*(N-2)^2/N/(N-1)}
    \cr
    Bootstrap
    \tab \eqn{S_P = S_0 + \sum_{i=1}^{S_0} (1 - p_i)^N}{S_P = S_0 + Sum
      (1-p_i)^N}
    }
    \code{specpool} normally uses basic Chao equation, but when there
    are no doubletons (\eqn{a2=0}) it switches to bias-corrected
    version. In that case the Chao equation simplifies to
    \eqn{S_0 + \frac{1}{2} a_1 (a_1-1) \frac{N-1}{N}}{S_0 + (N-1)/N * a1*(a1-1)/2}.

    The abundance-based estimates in \code{estimateR} use counts
    (numbers of individuals) of species in a single site. If called for
    a matrix or data frame, the function will give separate estimates
    for each site.  The two variants of extrapolated richness in
    \code{estimateR} are bias-corrected Chao and ACE (O'Hara 2005, Chiu
    et al. 2014).  The Chao estimate is similar as the bias corrected
    one above, but \eqn{a_i} refers to the number of species with
    abundance \eqn{i} instead of number of sites, and the small-sample
    correction is not used. The ACE estimate is defined as:

    \tabular{ll}{
    ACE
    \tab \eqn{S_P = S_{abund} + \frac{S_{rare}}{C_{ace}}+ \frac{a_1}{C_{ace}}
      \gamma^2_{ace}}{S_P = S_abund + S_rare/C_ace + a1/C_ace * gamma^2}
    \cr
    where \tab
    \eqn{C_{ace} = 1 - \frac{a_1}{N_{rare}}}{C_{ace} = 1- a1/N_{rare}}
    \cr
    \tab \eqn{\gamma^2_{ace} = \max \left[ \frac{S_{rare} \sum_{i=1}^{10}
      i(i-1)a_i}{C_{ace} N_{rare} (N_{rare} - 1)}-1, 0 \right]}{gamma^2 = 
    max(S_rare/C_ace (sum[i=1..10] i*(i-1)*a_i) / N_rare/(N_rare-1) -1 , 0)}
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
    no variance in the observed richness.  The equations of standard
    errors are too complicated to be reproduced in this help page, but
    they can be studied in the \R source code of the function and are
    discussed in the \code{\link{vignette}} \dQuote{diversity-vegan}
    that can be read with the \code{\link{vegandocs}} command. The
    standard error are based on the following sources: Chiu et
    al. (2014) for the Chao estimates and Smith and van Belle (1984)
    for the first-order Jackknife and the bootstrap (second-order
    jackknife is still missing).  For the variance estimator of
    \eqn{S_{ace}}{S_ace} see O'Hara (2005).

  Functions \code{poolaccum} and \code{estaccumR} are similar to
  \code{\link{specaccum}}, but estimate extrapolated richness indices
  of \code{specpool} or \code{estimateR} in addition to number of
  species for random ordering of sampling units. Function
  \code{specpool} uses presence data and \code{estaccumR} count
  data. The functions share \code{summary} and \code{plot}
  methods. The \code{summary} returns quantile envelopes of
  permutations corresponding the given level of \code{alpha} and
  standard deviation of permutations for each sample size. NB., these
  are not based on standard deviations estimated within \code{specpool}
  or \code{estimateR}, but they are based on permutations. The
  \code{plot} function shows the mean and envelope of permutations
  with given \code{alpha} for models. The selection of models can be
  restricted and order changes using the \code{display} argument in
  \code{summary} or \code{plot}. For configuration of \code{plot}
  command, see \code{\link[lattice]{xyplot}}.
}

\value{
  Function \code{specpool} returns a data frame with entries for
  observed richness and each of the indices for each class in
  \code{pool} vector.  The utility function \code{specpool2vect} maps
  the pooled values into a vector giving the value of selected
  \code{index} for each original site. Function \code{estimateR}
  returns the estimates and their standard errors for each
  site. Functions \code{poolaccum} and \code{estimateR} return
  matrices of permutation results for each richness estimator, the
  vector of sample sizes and a table of \code{means} of permutations
  for each estimator.
}

\references{
  Chao, A. (1987). Estimating the population size for capture-recapture
  data with unequal catchability. \emph{Biometrics} 43, 783--791.

  Chiu, C.H., Wang, Y.T., Walther, B.A. & Chao, A. (2014). Improved
  nonparametric lower bound of species richness via a modified
  Good-Turing frequency formula. \emph{Biometrics} 70, 671--682.
  
  Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial
  biodiversity through
  extrapolation. \emph{Phil. Trans. Roy. Soc. London} B 345, 101--118.

  O'Hara, R.B. (2005). Species richness estimators: how many species
  can dance on the head of a pin? \emph{J. Anim. Ecol.} 74, 375--386.

  Palmer, M.W. (1990). The estimation of species richness by
  extrapolation. \emph{Ecology} 71, 1195--1198.

  Smith, E.P & van Belle, G. (1984). Nonparametric estimation of
  species richness. \emph{Biometrics} 40, 119--129.
}
\author{Bob O'Hara (\code{estimateR}) and Jari Oksanen.}

\note{ The functions are based on assumption that there is a species
  pool: The community is closed so that there is a fixed pool size
  \eqn{S_P}.  In general, the functions give only the lower limit of
  species richness: the real richness is \eqn{S >= S_P}, and there is
  a consistent bias in the estimates. Even the bias-correction in Chao
  only reduces the bias, but does not remove it completely (Chiu et
  al. 2014).

  Optional small sample correction was added to \code{specpool} in
  \pkg{vegan} 2.2-0. It was not used in the older literature (Chao
  1987), but it is recommended recently (Chiu et al. 2014).

  See \url{http://viceroy.eeb.uconn.edu/EstimateS} for a more complete
  (and positive) discussion and alternative software for some platforms.
}
\seealso{\code{\link{veiledspec}}, \code{\link{diversity}}, \code{\link{beals}},
 \code{\link{specaccum}}. }
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
## Accumulation model
pool <- poolaccum(BCI)
summary(pool, display = "chao")
plot(pool)
## Quantitative model
estimateR(BCI[1:5,])
}
\keyword{ univar }

