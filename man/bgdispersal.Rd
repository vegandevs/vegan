\encoding{UTF-8}
\name{bgdispersal}
\alias{bgdispersal}

\title{ Coefficients of Biogeographical Dispersal Direction }

\description{ This function computes coefficients of dispersal direction
between geographically connected areas, as defined by Legendre and
Legendre (1984), and also described in Legendre and Legendre (1998,
section 13.3.4). }

\usage{
bgdispersal(mat, PAonly = FALSE, abc = FALSE)
}

\arguments{
  \item{mat}{ Data frame or matrix containing a community composition
    data table (species presence-absence or abundance data). }
  \item{PAonly}{ \code{FALSE} if the four types of coefficients, DD1 to
    DD4, are requested; \code{TRUE} if \code{DD1} and \code{DD2} only are
    sought (see Details). }
  \item{abc}{If \code{TRUE}, return tables \code{a}, \code{b} and \code{c}
    used in \code{DD1} and \code{DD2}.}
  }
\details{
The signs of the DD coefficients indicate the 
direction of dispersal, provided that the 
asymmetry is significant. A positive sign 
indicates dispersal from the first (row in DD 
tables) to the second region (column); a negative 
sign indicates the opposite. A McNemar test of 
asymmetry is computed from the presence-absence 
data to test the hypothesis of a significant 
asymmetry between the two areas under comparison.

In the input data table, the rows are sites or 
areas, the columns are taxa. Most often, the taxa 
are species, but the coefficients can be computed 
from genera or families as well. DD1 and DD2 only 
are computed for presence-absence data. The four 
types of coefficients are computed for 
quantitative data, which are converted to 
presence-absence for the computation of DD1 and 
DD2. \code{PAonly = FALSE} indicates that the four types 
of coefficients are requested. \code{PAonly = TRUE} if DD1 
and DD2 only are sought. }


\value{
Function \code{bgdispersal} returns a list containing the following matrices:
\item{ DD1 }{ \eqn{DD1[j,k] = (a * (b - c))/((a + b + c)^2)} }
\item{ DD2 }{ \eqn{DD2[j,k] = (2*a * (b - c))/((2*a + b + c) * (a + b +
    c))}
  where \eqn{a}, \eqn{b}, and \eqn{c} have the 
same meaning as in the computation of binary 
similarity coefficients. }
\item{ DD3 }{ DD3[j,k] = \eqn{W*(A-B) / ((A+B-W)^2)} }
\item{ DD4 }{ DD4[j,k] = \eqn{2*W*(A-B) / ((A+B)*(A+B-W)})
where \code{W = sum(pmin(vector1, vector2))}, \code{A = sum(vector1)},
\code{B = sum(vector2)} }
\item{ McNemar }{ McNemar chi-square statistic of asymmetry (Sokal and
Rohlf 1995): \eqn{2*(b*log(b) + c*log(c) - (b+c)*log((b+c)/2)) / q}
where \eqn{q = 1 + 1/(2*(b+c))} (Williams correction for continuity) }
\item{ prob.McNemar }{ probabilities associated 
with McNemar statistics, chi-square test. H0: no 
asymmetry in \eqn{(b-c)}. }
}


\references{ 
  Legendre, P. and V. Legendre. 1984. Postglacial dispersal of
  freshwater fishes in the Québec
  peninsula. \emph{Can. J. Fish. Aquat. Sci.} \strong{41}: 1781-1802.
  
  Legendre, P. and L. Legendre. 1998. \emph{Numerical ecology}, 2nd
  English edition. Elsevier Science BV, Amsterdam.
  
  Sokal, R. R. and F. J. Rohlf. 1995. \emph{Biometry. The principles and
  practice of statistics in biological research.} 3rd
  edn. W. H. Freeman, New York. 
}


\author{ Pierre Legendre, Departement de Sciences Biologiques,
  Universite de Montreal}

\note{The function uses a more powerful alternative for the McNemar test
  than the classical formula. The classical formula was constructed in
  the spirit of Pearson's Chi-square, but the formula in this function
  was constructed in the spirit of Wilks Chi-square or the \eqn{G}
  statistic. Function \code{\link{mcnemar.test}} uses the classical
  formula. The new formula was introduced in \pkg{vegan} version
  1.10-11, and the older implementations of \code{bgdispersal} used the
  classical formula.  }

\examples{
mat <- matrix(c(32,15,14,10,70,30,100,4,10,30,25,0,18,0,40,
  0,0,20,0,0,0,0,4,0,30,20,0,0,0,0,25,74,42,1,45,89,5,16,16,20),
  4, 10, byrow=TRUE)
bgdispersal(mat)
}

\keyword{ multivariate }
\keyword{ nonparametric }
