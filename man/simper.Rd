\encoding{UTF-8}
\name{simper}
\alias{simper}
\title{Similarity Percentages}

\description{
  Discriminating species between two groups using
  Bray-Curtis dissimilarities
}

\usage{
  simper(comm, group)
}
\arguments{
  \item{comm}{Community data}

  \item{group}{Vector assigning the groups}
}

\details{ Similarity percentage, \code{simper} (Clarke 1993) is based
  on the decomposition of Bray-Curtis dissimilarity index (see
  \code{\link{vegdist}}, \code{\link{designdist}}). The contribution
  of individual species \eqn{i} to the overall Bray-Curtis dissimilarity
  \eqn{d_{jk}}{d[jk]} is given by

  \deqn{d_{ijk} = \frac{|x_{ij}-x_{ik}|}{\sum_{i=1}^S (x_{ij}+x_{ik})}}{d[ijk] = abs(x[ij]-x[ik])/sum(x[ij]+x[ik])}
  
  where \eqn{x} is the abundance of species \eqn{i} in sampling units
  \eqn{j} and \eqn{k}. The overall index is the sum of the individual
  contributions over all \eqn{S} species 
  \eqn{d_{jk}=\sum_{i=1}^S d_{ijk}}{d[jk] = sum(i=1..S) d[ijk]}. 
  
  The \code{simper} functions performs pairwise comparisons of groups
  of sampling units and finds the average percentage contributions
  of each species to the average overall Bray-Curtis dissimilarity.
}

\value{
  A list of dataframes for every factor-combination.
  \item{contr}{average contribution to overall
  dissimilarity} \item{sd}{standard deviation of
  contribution} \item{contr/sd}{mean to sd ratio}
  \item{av_}{average abundance per group}
  \item{cum}{cumulative per cent contribution}
}

\examples{
data(dune)
data(dune.env)
with(dune.env, simper(dune, Management))
}
\author{
  Eduard Szöcs \email{szoe8822@uni-landau.de}
}

\references{
  Clarke, K.R. 1993. Non-parametric multivariate analyses of changes
    in community structure. \emph{Australian Journal of Ecology}, 18,
    117–143.
}
\keyword{multivariate}

