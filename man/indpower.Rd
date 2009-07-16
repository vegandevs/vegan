\encoding{UTF-8}
\name{indpower}
\alias{indpower}

\title{Indicator Power of Species}

\description{ Indicator power calculation of Halme et al.  (2009)
  scales up the congruence between indicator and target species.  


\usage{
indpower(x, type = 0)
}
\arguments{
  \item{x}{Community data frame or matrix. }
  \item{type}{Thehe type of
  statistic to be returned. See details for more explanation.}  
}
\details{

  Halme et al. (2009) described an index of indicator power defined as
  \code{IP_I = sqrt(a * b)}, where \code{a = S / O_I} and
  \code{a = 1 - (O_T - S) / (N - O_I)}. \code{N} is the number of sites,
  \code{S} is the number of shared occurrences of the indicator (\code{I})
  and the target (\code{T}) species. \code{O_I} and \code{O_T} are number
  of occurrences of the indicator and target species. The \code{type}
  argument in the function call enables to choose which statistic to
  return. \code{type = 0} returns \code{IP_I}, \code{type = 1} returns 
  \code{a}, \code{type = 2} returns \code{b}.
  Total indicator power (TIP) of an indicator species is the column mean
  (without its own value, see examples).  
  In the paper, one can find how to calculate confidence
  intervals for these statistics.

}
\value{
  A matrix with indicator species as rows and
  target species as columns (this is indicated by the first letters of the
  row/column names).
}
\references{
Halme, P., \enc{Mönkkönen}{Monkkonen}, M., Kotiaho, J. S, 
\enc{Ylisirniö}{Ylisirnio}, A-L. 2009. Quantifying the indicator power
of an indicator species. \emph{Conservation Biology} 23: 1008--1016.
}
\author{Peter Solymos}

\seealso{ \code{\link[labdsv]{duleg}} (package \pkg{labdsv}) for the
indicator species analysis of \enc{Dufrêne}{Dufrene} & Legendre. }

\examples{
data(dune)
## IP values
ip <- indpower(dune)
## and TIP values
diag(ip) <- NA
rowMeans(ip, na.rm=TRUE)
}
\keyword{ multivariate }

