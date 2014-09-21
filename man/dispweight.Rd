\encoding{UTF-8}
\name{dispweight}
\alias{dispweight}
\title{Dispersion-based weighting of species counts}

\description{Transform abundance data using dispersion-based weighting of species counts}

\usage{
dispweight(comm, groups, nsimul = 999, nullmodel = "c0_ind",
    plimit = 0.05)
}

\arguments{
  \item{comm}{Community data matrix.}
  \item{groups}{Factor describing the group structure.}
  \item{nsimul}{Number of simulations.}
  \item{nullmodel}{The \code{\link{nullmodel}} used in
    \code{\link{commsim}} within \code{groups}. The default
    follows Clarke et al. (2006).}
  \item{plimit}{Downweight species if their \eqn{p}-value is at or
     below this limit.}

}

\details{
Dispersion weighting can be used when there are \emph{a priori} defined groups 
(eg. lacation of samples). 

The dispersion index (\eqn{D}) is calculated as ratio between 
variance and mean and averaged across groups. D can then be used as divisor to 
downweight species abundances.

However species should be only downweighted, when there is evidence for 
over-dispersion (\eqn{D > 1}). This is tested using a permutational Chi-squared test. 
If this test shows statistically significant over-dispersion then the species is 
downweighted by \eqn{1/D}, else the divisor is \eqn{1}.
}

\value{
Function returns transformed data with the following new attributes:
    \item{D}{Average Dispersion index.}
    \item{p}{\eqn{p}-value of permutation test that \eqn{D = 1}.}
    \item{weights}{weights applied to community data.}
    \item{nsimul}{Number of simulations used to assess the \eqn{p}-value,
       or \code{NA} when simulations were not performed.}
}

\references{
Clarke, K. R., M. G. Chapman, P. J. Somerfield, and H. R. Needham. 2006. Dispersion-based Weighting of Species Counts in Assemblage Analyses. \emph{Marine Ecology Progress Series}, 320, 11–27.
}
\author{
  Eduard Szöcs \email{szoe8822@uni-landau.de}
}


\examples{
data(dune)
data(dune.env)
## transformed community data
dpw <- dispweight(dune, dune.env$Management, nsimul = 99)
}


\keyword{multivariate}

