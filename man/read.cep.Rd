\name{read.cep}
\alias{read.cep}

\title{Reads a CEP (Canoco) data file }
\description{
  \code{read.cep} reads a file formatted by relaxed strict CEP format
  used by \code{Canoco} software, among others.
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
  designed for punched cards.  One of these was the `condensed strict'
  format which was adopted by popular software \code{DECORANA} and
  \code{TWINSPAN}.

  The condensed CEP and CANOCO formats have:
  \itemize{
    \item Two or three title cards, most importantly specifying the format
    and the number of items per record.
    \item Data in condensed format: First number on the line is the site
      identifier, and it is followed by pairs (`couplets') of numbers
      identifying the species and its abundance (an integer and a floating
      point number).
    \item Species and site names, given in Fortran format \code{(10A8)}:
      Ten names per line, eight columns for each.
  }

  With option \code{positive = TRUE} the function removes all lines and
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
  Ter Braak, C.J.F. (1984--): CANOCO -- a FORTRAN program for \emph{cano}nical
  \emph{c}ommunity \emph{o}rdination by [partial] [detrended] [canonical]
  correspondence analysis, principal components analysis and redundancy
  analysis. \emph{TNO Inst. of Applied Computer Sci., Stat. Dept. Wageningen,
  The Netherlands}. 
}
\author{ Jari Oksanen }

\note{

  Fortran function was used to read data in \pkg{vegan} 2.4-5 and
  earlier, but Fortran I/O is no longer allowed in CRAN packages, and
  this function was re-written in \R. The original Fortran code was
  more robust, and there are several legacy data sets that may fail
  with the current version, but could be read with the previous
  Fortran version. Github package \pkg{cepreader}
  \url{https://github.com/vegandevs/cepreader/} makes available the
  original Fortran code. The \pkg{cepreader} package can also read
  \sQuote{free} and \sQuote{open} CANOCO formats that are not handled
  with this function.

}

\section{Warning}{

  The function is based on \code{\link[utils]{read.fortran}} which can
  do really bad job in interpreting Fortran format, in particular,
  when abundance data were defined with decimal format.  For instance,
  when format \code{F5.1} specifies one decimal for abundances,
  function \code{read.fortran} will divide data by ten, and input 0.5
  becomes 0.05 and 50 becomes 5. You should check the scaling of
  results after reading the data, and if necessary, multiply results
  to the original scale.

}


\examples{
## Provided that you have the file `dune.spe'
\dontrun{
theclassic <- read.cep("dune.spe", force=T)}
}
\keyword{ IO }
\keyword{ file }




