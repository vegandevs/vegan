\name{make.cepnames}
\alias{make.cepnames}

\title{Abbreviates a Two-Part Botanical or Zoological Latin Name into
  Character String}

\description{ Function is based on \code{\link{abbreviate}}, and will
  take given number of characters from the first (genus) and last
  (epithet) component of botanical or zoological Latin name and combine
  these into one shorter character string. The names will be unique and
  more characters will be used if needed. The default usage makes names
  with 4+4 characters popularized in Cornell Ecology Programs (CEP) and
  often known as CEP names. Abbreviated names are useful in ordination
  plots and other graphics to reduce clutter.  }

\usage{
make.cepnames(names, minlengths = c(4,4), seconditem = FALSE,
   uniqgenera = FALSE, named = FALSE, method)
}

\arguments{
  \item{names}{The names to be abbreviated into a vector abbreviated names.}
  \item{minlengths}{The minimum lengths of first and second part of the
    abbreviation. If abbreviations are not unique, the parts can be
    longer.}
  \item{seconditem}{Take always the second part of the original name
    to the abbreviated name instead of the last part.}
  \item{uniqgenera}{Should the first part of the abbreviation (genus)
    also be unique. Unique genus can take space from the second part
    (epithet).}
  \item{method}{The \code{\link{abbreviate}} argument in last attempt to
    abbreviate the abbreviation. The default \code{method} tries
    to drop character from the end, but \code{"both.sides"} can remove
    characters from any position, including the genus part, and same
    genus can be abbreviated differently.}
  \item{named}{Should the result vector be named by original
    \code{names}.}
}

\details{ Cornell Ecology Programs (CEP) used eight-letter abbreviations
  for species and site names. In species, the names were formed by
  taking four first letters of the generic name and four first letters
  of the specific or subspecific epithet. The current function produces
  4+4 character CEP names as default, but it can also use other
  lengths. The function is based on \code{\link{abbreviate}} and will
  produce longer names if constructed names are not unique. If generic
  name is shorter than specified minimum length, more characters can be
  used by the epithet. If \code{uniqgenera = TRUE}, unique generic name
  may need more characters than the minimum length, and these reduce the
  number of characters available for the epithet. The function drops
  characters from the end, but with \code{method = "both.sides"} the
  function tries to drop characters from other positions in the final
  attempt to abbreviate abbreviations, starting with lower-case
  vowels. The constructed names are based on binomial names with two
  components. If there is only one component (typically a genus), it can
  use also the space reserved for the second component. If there are
  more than two components of input names, the function selects only
  two: always the first, and either the last or second depending on the
  argument \code{seconditem}. Function will issue a warning of
  duplicated input names or of duplicated inferred binomial names, but
  resulting CEP names are made unique by adding a sequence number (with
  a warning listing these names). }

\value{
  Function returns a vector of abbreviated names.
}
\author{ Jari Oksanen }

\seealso{
 \code{\link{abbreviate}}.
}

\note{The function does not handle Author names except strictly
  two-part names with \code{seconditem = TRUE}. It is often necessary to
  edit abbreviations manually.}

\examples{
names <- c("Aa maderoi", "Capsella bursa-pastoris", "Taraxacum",
  "Cladina rangiferina", "Cladonia rangiformis", "Cladonia cornuta",
  "Cladonia cornuta var. groenlandica", "Euphorbia cf. enopla",
  "Rumex acetosa", "Rumex acetosella")
make.cepnames(names)
make.cepnames(names, uniqgenera = TRUE)
make.cepnames(names, method = "both.sides")
}
\keyword{ character }

