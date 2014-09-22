\encoding{UTF-8}
\name{dispweight}
\alias{dispweight}
\alias{gdispweight}
\title{Dispersion-based weighting of species counts}

\description{Transform abundance data using dispersion-based weighting
of species counts}

\usage{
dispweight(comm, groups, nsimul = 999, nullmodel = "c0_ind",
    plimit = 0.05)
gdispweight(formula, data, plimit = 0.05)
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
  \item{formula, data}{Formula where the left-hand side is the
    community data frame and right-hand side gives the explanatory
    variables. The explanatory variables are found in the data frame
    given in \code{data}.}

}

\details{

Dispersion weighting can be used when there are \emph{a priori}
defined groups (eg. location of samples). The dispersion index
(\eqn{D}) is calculated as ratio between variance and mean and
averaged across groups. \eqn{D} can then be used as divisor to
downweight species abundances.  However species should be only
downweighted, when there is evidence for over-dispersion (\eqn{D >
1}). This is tested using a permutational Chi-squared test.  If this
test shows statistically significant over-dispersion then the species
is downweighted by \eqn{1/D}, else the divisor is \eqn{1} (Clarke et
al. 2006).

Function \code{gdispweight} is a generalized parametric version
dispersion weighting of Clarke et al. (2006). The original method
analysed Pearson Chi-square with respect to a single factor, but the
generalized function is based on \code{\link{glm}} with
\code{\link{quasipoisson}} error \{code{\link{family}}. Any
\code{\link{glm}} model can be used, including continuous
covariates. The original method estimates significance using
multinomial simulation for each species by groups
(\code{link{commsim}} model \code{"c0_ind"}), but the generalized
function uses parametric test of the Pearson Chi-square.  Further, the
residual degrees of freedom are estimated differently: the original
method ignores factor levels where species is absent, but the
generalized function estimates the degrees of freedom as in
\code{\link{glm}}. Apart from different way of estimating degrees of
freedom, the downweighting is identical to the original method.

}

\value{
Function returns transformed data with the following new attributes:
    \item{D}{Dispersion index.}
    \item{df}{Degrees of freedom for each species.}
    \item{p}{\eqn{p}-value of permutation test that \eqn{D = 1}.}
    \item{weights}{weights applied to community data.}
    \item{nsimul}{Number of simulations used to assess the \eqn{p}-value,
       or \code{NA} when simulations were not performed.}
}

\references{
Clarke, K. R., M. G. Chapman, P. J. Somerfield, and
H. R. Needham. 2006. Dispersion-based weighting of species counts in
assemblage analyses. \emph{Marine Ecology Progress Series}, 320,
11–27.
}
\author{
  Eduard Szöcs \email{eduardszoesc@gmail.com} and Jari Oksanen.
}


\examples{
data(dune)
data(dune.env)
## transformed community data
dunepw <- with(dune.env, dispweight(dune, Management, nsimul = 99))
dunegpw <- gdispweight(dune ~ Management, dune.env)
}


\keyword{multivariate}

