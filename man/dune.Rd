\name{dune}
\alias{dune}
\alias{dune.env}
\docType{data}
\title{Vegetation and Environment in Dutch Dune Meadows. }
\usage{
  data(dune)
  data(dune.env)
}
\description{
  The dune meadow vegetation data, \code{dune}, has cover class values
  of 30 species on 20 sites. The corresponding environmental data frame
  \code{dune.env} has following entries: } 

\format{
  \code{dune} is a data frame of observations of 30 species at 20
  sites. The species names are abbreviated to 4+4 letters (see
  \code{\link{make.cepnames}}). The following names are changed from
  the original source (Jongman et al. 1987): \emph{Leontodon
  autumnalis} to \emph{Scorzoneroides}, and \emph{Potentilla
  palustris} to \emph{Comarum}.
  
  \code{dune.env} is a data frame of 20 observations on the following
  5 variables:
  \describe{
    \item{A1:}{a numeric vector of thickness of soil A1 horizon.}
    \item{Moisture:}{an ordered factor with levels: \code{1} < \code{2} <
      \code{4} < \code{5}.}
    %\item{Moisture}{\code{1} < \code{2} < \code{4} < \code{5}}
    \item{Management:}{a factor with levels: \code{BF} (Biological
      farming), \code{HF} (Hobby farming), \code{NM} (Nature
      Conservation Management), and \code{SF} (Standard Farming).}
    %\item{Management}{\code{BF}: Biological Farming  }
    %\item{Management}{\code{HF}: Hobby Farming }
    %\item{Management}{\code{NM}: Nature Conservation Management }
    %\item{Management}{\code{SF}: Standard Farming }
    \item{Use:}{an ordered factor of land-use with levels: \code{Hayfield}
      < \code{Haypastu} < \code{Pasture}.}
    %\item{Use}{\code{Hayfield} < \code{Haypastu} < \code{Pasture}}
    \item{Manure:}{an ordered factor with levels: \code{0} < \code{1} <
      \code{2} < \code{3} < \code{4}.} 
    %\item{Manure}{\code{0} < \code{1} < \code{2} < \code{3} < \code{4}}
  }
}
\source{
  Jongman, R.H.G, ter Braak, C.J.F & van Tongeren,
  O.F.R. (1987). \emph{Data Analysis in Community and Landscape
    Ecology}. Pudoc, Wageningen.
}
\examples{
data(dune)
data(dune.env)
}
\keyword{datasets}
