\encoding{UTF-8}
\name{BCI}
\alias{BCI}
\docType{data}
\title{Barro Colorado Island Tree Counts}
\description{
  Tree counts in 1-hectare plots in the Barro Colorado Island.
}
\usage{data(BCI)}
\format{
  A data frame with 50 plots (rows) of 1 hectare with counts of trees on each
  plot with total of 225 species (columns). Full Latin names are used
  for tree species.
 }
\details{
  Data give the numbers of trees at least 10 cm in
  diameter at breast height (1.3 m above the ground) in each one hectare
  square of forest. Within each one hectare square, all individuals of
  all species were tallied and are recorded in this table.

  The data frame contains only the Barro Colorado Island subset of the
  original data.

  The quadrats are located in a regular grid. See\code{examples} for the
  coordinates. 

}
\source{
  \url{http://www.sciencemag.org/cgi/content/full/295/5555/666/DC1}
}
\references{
  Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster,
  R.B., \enc{Nuñez}{Nunez}, P., Aguilar, S., Valencia, R., Villa, G., Muller-Landau,
  H.C., Losos, E. & Hubbell, S.P. (2002). Beta-diversity in tropical
  forest trees. \emph{Science} 295, 666--669.
}
\seealso{\code{\link[BiodiversityR]{BCI.env}} in \pkg{BiodiversityR}
  package for environmental data  (coordinates are given below in the
  \code{examples}).}  
\examples{
data(BCI)
## UTM Coordinates (in metres)
UTM.EW <- rep(seq(625754, 626654, by=100), each=5)
UTM.NS <- rep(seq(1011569,  1011969, by=100), len=50)
}
\keyword{datasets}
