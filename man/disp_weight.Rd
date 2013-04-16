\encoding{UTF-8}
\name{disp_weight}
\alias{disp_weight}
\title{Dispersion-based weighting of species counts}

\description{Transform abundance data using dispersion-based weighting of species counts}

\usage{disp_weight(comm, group, nperm = 1000)}

\arguments{
  \item{comm}{Community data matrix.}
  \item{group}{Factor describing the group structure.}
  \item{nperm}{Number of permutations}
}

\details{
Dispersion weighting can be used when there are \emph{a priori} defined groups 
(eg. lacation of samples). 

The dispersion index (D) is calculated as ratio between 
variance and mean and averaged across groups. D can then be used as divisor to 
downweight species abundances.

However species should be only downweighted, when there is evidence for 
over-dispersion (D > 1). This is tested using a permutational Chi-squared test. 
If this test shows statistically significant over-dispersion then the species is 
downweighted by 1/D, else the divisor is 1.
}

\value{
A data.frame with three columns:
    \item{D}{Average Dispersion index}
    \item{p}{p-value of permutation test that D = 1}
    \item{weights}{weights to apply to community data}
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
# calculate weights
dpw <- disp_weight(dune, dune.env$Management, nperm = 100)
# calculate transformed data
sweep(dune, MARGIN = 2, dpw$weights, `*`)
}


\keyword{multivariate}

