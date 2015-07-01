\name{sipoo}
\alias{sipoo}
\docType{data}
\title{ Birds in the Archipelago of Sipoo (Sibbo)}
\description{
Land birds on islands covered by 
coniferous forest in the Sipoo archipelago, southern Finland (landbridge/
oceanic distinction unclear from source). 
}
\usage{data(sipoo)}
\format{
  A data frame with 18 sites and 50 species (Simberloff & Martin, 1991,
  Appendix 3). The species are referred by 4+4 letter abbreviation of
  their Latin names (but using five letters in two species names to make
  these unique). The example gives the areas of the studies islands in
  hectares. 
}

\source{
Simberloff, D. & Martin, J.-L.  (1991). 
Nestedness of insular avifaunas:
simple summary statistics masking
complex species patterns.  \emph{Ornis
Fennica} 68:178--192. 
}
\examples{
data(sipoo)
## Areas of the islands in hectares
sipoo.area <-  c(1.1, 2.1, 2.2, 3.1, 3.5, 5.8, 6, 6.1, 6.5, 11.4, 13,
14.5, 16.1 ,17.5, 28.7, 40.5, 104.5, 233) 
}
\keyword{datasets}
