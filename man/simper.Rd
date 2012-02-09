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
  Clarke, K.R. 1993. Non-parametric multivariate analyses
    of changes in community structure. \emph{Austral Ecology},
    18, 117–143.
}
\keyword{multivariate}

