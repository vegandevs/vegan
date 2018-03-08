\encoding{UTF-8}
\name{nestedtemp}
\alias{nestedtemp}
\alias{nestedchecker}
\alias{nestedn0}
\alias{nesteddisc}
\alias{nestednodf}
\alias{nestedbetasor}
\alias{nestedbetajac}
\alias{plot.nestedtemp}
\alias{plot.nestednodf}

\title{ Nestedness Indices for Communities of Islands or Patches }

\description{
  Patches or local communities are regarded as nested if they all could
  be subsets of the same community. In general, species poor communities
  should be subsets of species rich communities, and rare species should
  only occur in species rich communities.
}

\usage{
nestedchecker(comm)
nestedn0(comm)
nesteddisc(comm, niter = 200)
nestedtemp(comm, ...)
nestednodf(comm, order = TRUE, weighted = FALSE, wbinary = FALSE)
nestedbetasor(comm)
nestedbetajac(comm)
\method{plot}{nestedtemp}(x, kind = c("temperature", "incidence"),
    col=rev(heat.colors(100)),  names = FALSE, ...)
\method{plot}{nestednodf}(x, col = "red", names = FALSE, ...)
}

\arguments{
  \item{comm}{Community data.}
  \item{niter}{Number of iterations to reorder tied columns.}
  \item{x}{Result object for a \code{plot}.}
  \item{col}{Colour scheme for matrix temperatures.}
  \item{kind}{The kind of plot produced.}
  \item{names}{Label columns and rows in the plot using names in \code{comm}.
    If it is a logical vector of length 2, row and column labels are
    returned accordingly.}
  \item{order}{Order rows and columns by frequencies.}
  \item{weighted}{Use species abundances as weights of interactions.}
  \item{wbinary}{Modify original method so that binary data give the same
    result in weighted and and unweighted analysis. }
  \item{\dots}{Other arguments to functions.}
}

\details{
  The nestedness functions evaluate alternative indices of nestedness.
  The functions are intended to be used together with Null model
  communities and used as an argument in \code{\link{oecosimu}} to analyse
  the non-randomness of results.

  Function \code{nestedchecker} gives the number of checkerboard units,
  or 2x2 submatrices where both species occur once but on different
  sites (Stone & Roberts 1990).

  Function \code{nestedn0} implements
  nestedness measure N0 which is the number of absences from the sites
  which are richer than the most pauperate site species occurs
  (Patterson & Atmar 1986).

  Function \code{nesteddisc} implements discrepancy index which is the
  number of ones that should be shifted to fill a row with ones in a
  table arranged by species frequencies (Brualdi & Sanderson
  1999). The original definition arranges species (columns) by their
  frequencies, but did not have any method of handling tied
  frequencies.  The \code{nesteddisc} function tries to order tied
  columns to minimize the discrepancy statistic but this is rather
  slow, and with a large number of tied columns there is no guarantee
  that the best ordering was found (argument \code{niter} gives the
  maximum number of tried orders). In that case a warning of tied
  columns will be issued.

  Function \code{nestedtemp} finds the matrix temperature which is
  defined as the sum of \dQuote{surprises} in arranged matrix.  In
  arranged unsurprising matrix all species within proportion given by
  matrix fill are in the upper left corner of the matrix, and the
  surprise of the absence or presences is the diagonal distance from the
  fill line (Atmar & Patterson 1993). Function tries to pack species and
  sites to a low temperature (\enc{Rodríguez-Gironés}{Rodriguez-Girones}
  & Santamaria 2006), but this is an iterative procedure, and the
  temperatures usually vary among runs.  Function \code{nestedtemp} also
  has a \code{plot} method which can display either incidences or
  temperatures of the surprises. Matrix temperature was rather vaguely
  described (Atmar & Patterson 1993), but
  \enc{Rodríguez-Gironés}{Rodriguez-Girones} & Santamaria (2006) are
  more explicit and their description is used here. However, the results
  probably differ from other implementations, and users should be
  cautious in interpreting the results. The details of calculations are
  explained in the \code{\link{vignette}} \emph{Design decisions and
  implementation} that you can read using functions
  \code{\link{browseVignettes}}. Function
  \code{\link[bipartite]{nestedness}} in the \pkg{bipartite} package is
  a direct port of the BINMATNEST programme of
  \enc{Rodríguez-Gironés}{Rodriguez-Girones} & Santamaria (2006).

  Function \code{nestednodf} implements a nestedness metric based on
  overlap and decreasing fill (Almeida-Neto et al., 2008). Two basic
  properties are required for a matrix to have the maximum degree of
  nestedness according to this metric: (1) complete overlap of 1's
  from right to left columns and from down to up rows, and (2)
  decreasing marginal totals between all pairs of columns and all
  pairs of rows. The nestedness statistic is evaluated separately for
  columns (\code{N columns}) for rows (\code{N rows}) and combined for
  the whole matrix (\code{NODF}).  If you set \code{order = FALSE},
  the statistic is evaluated with the current matrix ordering allowing
  tests of other meaningful hypothesis of matrix structure than
  default ordering by row and column totals (breaking ties by total
  abundances when \code{weighted = TRUE}) (see Almeida-Neto et
  al. 2008). With \code{weighted = TRUE}, the function finds the
  weighted version of the index (Almeida-Neto & Ulrich,
  2011). However, this requires quantitative null models for adequate
  testing. Almeida-Neto & Ulrich (2011) say that you have positive
  nestedness if values in the first row/column are higher than in the
  second.  With this condition, weighted analysis of binary data will
  always give zero nestedness. With argument \code{wbinary = TRUE},
  equality of rows/columns also indicates nestedness, and binary data
  will give identical results in weighted and unweighted analysis.
  However, this can also influence the results of weighted analysis so
  that the results may differ from Almeida-Neto & Ulrich (2011).

  Functions \code{nestedbetasor} and \code{nestedbetajac} find
  multiple-site dissimilarities and decompose these into components of
  turnover and nestedness following Baselga (2012); the pairwise
  dissimilarities can be found with \code{\link{designdist}}. This can
  be seen as a decomposition of beta diversity (see
  \code{\link{betadiver}}).  Function \code{nestedbetasor} uses
  \enc{Sørensen}{Sorensen} dissimilarity and the turnover component is
  Simpson dissimilarity (Baselga 2012), and \code{nestedbetajac} uses
  analogous methods with the Jaccard index. The functions return a
  vector of three items: turnover, nestedness and their sum which is
  the multiple \enc{Sørensen}{Sorensen} or Jaccard dissimilarity. The
  last one is the total beta diversity (Baselga 2012). The functions
  will treat data as presence/absence (binary) and they can be used
  with binary \code{\link{nullmodel}}. The overall dissimilarity is
  constant in all \code{\link{nullmodel}}s that fix species (column)
  frequencies (\code{"c0"}), and all components are constant if row
  columns are also fixed (e.g., model \code{"quasiswap"}), and the
  functions are not meaningful with these null models.

}

\value{
  The result returned by a nestedness function contains an item called
  \code{statistic}, but the other components differ among functions. The
  functions are constructed so that they can be handled by
  \code{\link{oecosimu}}.
}

\references{

  Almeida-Neto, M., \enc{Gumarães}{Gumaraes}, P.,
  \enc{Gumarães}{Gumaraes}, P.R., Loyola, R.D. & Ulrich, W. (2008). A
  consistent metric for nestedness analysis in ecological systems:
  reconciling concept and measurement. \emph{Oikos} 117, 1227--1239.

  Almeida-Neto, M. & Ulrich, W. (2011). A straightforward
  computational approach for measuring nestedness using quantitative
  matrices. \emph{Env. Mod. Software} 26, 173--178.

  Atmar, W. & Patterson, B.D. (1993). The measurement of order and
  disorder in the distribution of species in fragmented
  habitat. \emph{Oecologia} 96, 373--382.

  Baselga, A. (2012). The relationship between species replacement,
  dissimilarity derived from nestedness, and nestedness. \emph{Global
  Ecol. Biogeogr.} 21, 1223--1232.

  Brualdi, R.A. & Sanderson, J.G. (1999). Nested species subsets, gaps,
  and discrepancy. \emph{Oecologia} 119, 256--264.

  Patterson, B.D. & Atmar, W. (1986). Nested subsets and the structure
  of insular mammalian faunas and archipelagos. \emph{Biol. J. Linnean
    Soc.} 28, 65--82.

  \enc{Rodríguez-Gironés}{Rodriguez-Girones}, M.A.  & Santamaria, L.
  (2006). A new algorithm to calculate the nestedness temperature of
  presence-absence matrices. \emph{J. Biogeogr.} 33, 924--935.

  Stone, L. & Roberts, A. (1990). The checkerboard score and species
  distributions. \emph{Oecologia} 85, 74--79.

  Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar,
  W. (1998). A comparative analysis of nested subset patterns of species
  composition. \emph{Oecologia} 113, 1--20.
  }
\author{ Jari Oksanen and Gustavo Carvalho (\code{nestednodf}). }

\seealso{
  In general, the functions should be used with \code{\link{oecosimu}}
  which generates Null model communities to assess the non-randomness of
  nestedness patterns.
 }
\examples{
data(sipoo)
## Matrix temperature
out <- nestedtemp(sipoo)
out
plot(out)
plot(out, kind="incid")
## Use oecosimu to assess the non-randomness of checker board units
nestedchecker(sipoo)
oecosimu(sipoo, nestedchecker, "quasiswap")
## Another Null model and standardized checkerboard score
oecosimu(sipoo, nestedchecker, "r00", statistic = "C.score")
}
\keyword{ univar }

