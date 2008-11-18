\encoding{UTF-8}
\name{varespec}
\alias{varechem}
\alias{varespec}
\docType{data}
\title{Vegetation and environment in lichen pastures}
\usage{
       data(varechem)
       data(varespec)
}
\description{
  The \code{varespec} data frame has 24 rows and 44 columns.  Columns
  are estimated cover values of 44 species.  The variable names are
  formed from the scientific names, and are self explanatory for anybody
  familiar with the vegetation type.
The \code{varechem} data frame has 24 rows and 14 columns, giving the
soil characteristics of the very same sites as in the \code{varespec}
data frame. The chemical measurements have obvious names.
\code{Baresoil} gives the estimated cover of bare soil, \code{Humpdepth}
the thickness of the humus layer.
}
\references{
\enc{Väre}{Vare}, H., Ohtonen, R. and Oksanen, J. (1995) Effects of reindeer
grazing on understorey vegetation in dry Pinus sylvestris
forests. \emph{Journal of Vegetation Science} 6, 523--530.  
}
\examples{
data(varespec)
data(varechem)
}
\keyword{datasets}
