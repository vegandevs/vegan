\name{make.cepnames}
\alias{make.cepnames}

\title{Abbreviates a Botanical or Zoological Latin Name into an Eight-character Name} 
\description{
  A standard CEP name has four first letters of the generic name and
  four first letters of the specific epithet of a Latin name. The last
  epithet, that may be a subspecific name, is used in the current
  function. If the name has only one component, it is abbreviated to
  eight characters (see \code{\link{abbreviate}}).
  The returned names are made unique with function
  \code{\link{make.unique}} which adds numbers to the end of CEP names if needed.
}
\usage{
make.cepnames(names)
}
\arguments{
  \item{names}{The names to be formatted into CEP names. }
}
\details{
  Cornell Ecology Programs (CEP) used eight-letter abbreviations for
  species and site names. In species, the names were formed by taking
  four first letters of the generic name and four first letters of the
  specific or subspecific epithet. The CEP names were originally used,
  because old \code{FORTRAN IV} did not have \code{CHARACTER} data type,
  but text  had to be stored in numerical variables, which in
  popular computers could hold four characters. In modern times,
  there is no reason for this limitation, but ecologists are used to
  these names, and they may be practical to avoid congestion in
  ordination plots. 
}
\value{
  Function returns CEP names.
}
\author{ Jari Oksanen }
\note{
  The function is simpleminded and rigid. You must write a better one if
  you need.  
}
\seealso{
  \code{\link{make.names}}, \code{\link{strsplit}},
  \code{\link{substring}}, \code{\link{paste}}, \code{\link{abbreviate}}.
  }
\examples{
make.cepnames(c("Aa maderoi", "Poa sp.", "Cladina rangiferina",
"Cladonia cornuta", "Cladonia cornuta var. groenlandica",
"Cladonia rangiformis", "Bryoerythrophyllum"))
data(BCI)
colnames(BCI) <- make.cepnames(colnames(BCI))
}
\keyword{ character }

