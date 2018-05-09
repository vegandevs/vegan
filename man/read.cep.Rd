\name{read.cep}
\alias{read.cep}

\title{Reads a CEP (Canoco) data file }
\description{
  \code{read.cep} reads a file formatted with relaxed strict CEP format
  used in Canoco software, among others.
}
\usage{
read.cep(file, positive=TRUE)
}

\arguments{
  \item{file}{File name (character variable). }
  \item{positive}{Only positive entries, like in community data.}
}
\details{

  Cornell Ecology Programs (CEP) introduced several data formats
  designed for punched cards.  One of these was the \sQuote{condensed
  strict} format which was adopted by popular software DECORANA and
  TWINSPAN. A relaxed variant of this format was later adopted in
  Canoco software (ter Braak 1984). Function \code{read.cep} reads
  legacy files written in this format.

  The condensed CEP and CANOCO formats have:
  \itemize{
    \item Two or three title cards, most importantly specifying the format
      and the number of items per record.
    \item Data in condensed format: First number on the line is the
      site identifier (an integer), and it is followed by pairs
      (\sQuote{couplets}) of numbers identifying the species and its
      abundance (an integer and a floating point number).
    \item Species and site names, given in Fortran format \code{(10A8)}:
      Ten names per line, eight columns for each.
  }

  With option \code{positive = TRUE} the function removes all rows and
  columns with zero or negative marginal sums.  In community data
  with only positive entries, this removes empty sites and species.
  If data entries can be negative, this ruins data, and such data sets
  should be read in with option \code{positive = FALSE}.
}
\value{
  Returns a data frame, where columns are species and rows are
  sites. Column and row names are taken from the CEP file, and changed
  into unique \R names by \code{\link{make.names}} after stripping the blanks.
}
\references{ 
  ter Braak, C.J.F. (1984--): CANOCO -- a FORTRAN program for \emph{cano}nical
  \emph{c}ommunity \emph{o}rdination by [partial] [detrended] [canonical]
  correspondence analysis, principal components analysis and redundancy
  analysis. \emph{TNO Inst. of Applied Computer Sci., Stat. Dept. Wageningen,
  The Netherlands}. 
}
\author{ Jari Oksanen }

\note{

  Function \code{read.cep} used Fortran to read data in \pkg{vegan}
  2.4-5 and earlier, but Fortran I/O is no longer allowed in CRAN
  packages, and the function was re-written in \R. The original
  Fortran code was more robust, and there are several legacy data sets
  that may fail with the current version, but could be read with the
  previous Fortran version. CRAN package \pkg{cepreader} makes
  available the original Fortran-based code run in a separate
  subprocess. The \pkg{cepreader} package can also read \sQuote{free}
  and \sQuote{open} Canoco formats that are not handled in this
  function.

  The function is based on \code{\link[utils]{read.fortran}}. If the
  \code{REAL} format defines a decimal part for species abundances
  (such as \code{F5.1}), \code{\link[utils]{read.fortran}} divides the
  input with the corresponding power of 10 even when the input data
  had explicit decimal separator. With \code{F5.1}, 100 would become
  10, and 0.1 become 0.01. Function \code{read.cep} tries to undo this
  division, but you should check the scaling of results after reading
  the data, and if necessary, multiply results to the original scale.

}


\examples{
## Provided that you have the file "dune.spe"
\dontrun{
theclassic <- read.cep("dune.spe")}
}
\keyword{ IO }
\keyword{ file }




