\name{mite}
\alias{mite}
\alias{mite.env}
\alias{mite.pcnm}
\alias{mite.xy}
\docType{data}
\title{Oribatid Mite Data with Explanatory Variables
  }
\description{
Oribatid mite data. 70 soil cores collected by Daniel Borcard in 1989.
See Borcard et al. (1992, 1994) for details.
}
\usage{
data(mite)
data(mite.env)
data(mite.pcnm)
data(mite.xy)
}
\format{
  There are three linked data sets: \code{mite} that contains the data
  on 35 species of Oribatid mites, \code{mite.env} that contains
  environmental data in the same sampling sites, \code{mite.xy}
  that contains geographic coordinates,  and \code{mite.pcnm}
  that contains 22 PCNM base functions (columns) computed from the geographic
  coordinates of the 70 sampling sites (Borcard & Legendre 2002).
  The whole sampling area was 2.5 m x 10 m in size.

  The fields in the environmental data are:
  \describe{
    \item{SubsDens}{Substrate density (g/L)}
    \item{WatrCont}{Water content of the substrate (g/L)}
    \item{Substrate}{Substrate type, factor with levels \code{Sphagn1,
	Sphagn2 Sphagn3 Sphagn Litter Barepeat Interface} }
    \item{Shrub}{Shrub density, an ordered factor with levels \code{1} <
      \code{2} < \code{3}}
    \item{Topo}{Microtopograhy, a factor with levels \code{Blanket} and \code{Hummock}}
  }
}
\source{
  Pierre Legendre
}
\references{
Borcard, D., P. Legendre and P. Drapeau. 1992. Partialling out the
spatial component of ecological variation. Ecology 73: 1045-1055.

Borcard, D. and P. Legendre. 1994. Environmental control and spatial
structure in ecological communities: an example using Oribatid mites
(Acari, Oribatei). Environmental and Ecological Statistics 1: 37-61.

Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of
ecological data by means of principal coordinates of neighbour
matrices. Ecological Modelling 153: 51-68.
}
\examples{
data(mite)
}
\keyword{datasets}
