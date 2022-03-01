\encoding{UTF-8}
\name{sipoo}
\alias{sipoo}
\alias{sipoo.map}
\docType{data}
\title{ Birds in the Archipelago of Sipoo (Sibbo and Borgå)}
\description{
Land birds on islands covered by 
coniferous forest in the Sipoo Archipelago, southern Finland.
}
\usage{
  data(sipoo)
  data(sipoo.map)
}

\format{
  The \code{sipoo} data frame contains data of occurrences of 50 land
  bird species on 18 islands in the Sipoo Archipelago (Simberloff &
  Martin, 1991, Appendix 3). The species are referred by 4+4 letter
  abbreviation of their Latin names (but using five letters in two
  species names to make these unique).

  The \code{sipoo.map} data contains the geographic coordinates of the
  islands in the ETRS89-TM35FIN coordinate system (EPSG:3067) and the
  areas of islands in hectares.
}

\source{
  Simberloff, D. & Martin, J.-L.  (1991).  Nestedness of insular
  avifaunas: simple summary statistics masking complex species patterns.
  \emph{Ornis Fennica} 68:178--192.
}

\examples{
data(sipoo)
data(sipoo.map)
plot(N ~ E, data=sipoo.map, asp = 1)
}

\keyword{datasets}
