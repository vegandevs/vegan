\name{mantel.correlog}
\alias{mantel.correlog}
\alias{print.mantel.correlog}
\alias{plot.mantel.correlog}
\title{ Mantel correlogram }
\description{
Function \code{mantel.correlog} computes a multivariate Mantel correlogram. 
}
\usage{
mantel.correlog(D.eco, D.geo=NULL, XY=NULL, n.class=0, break.pts=NULL, 
cutoff=TRUE, r.type="pearson", nperm=999)

\method{plot}{mantel.correlog}(x,...)
}

\arguments{
  \item{D.eco}{ An ecological distance matrix, with class either \code{dist} or \code{matrix}. }
  \item{D.geo}{ A geographic distance matrix, with class either \code{dist} or \code{matrix}. Provide either \code{D.geo} or \code{XY}. Default: \code{D.geo=NULL}. }
  \item{XY}{ A file of Cartesian geographic coordinates of the points. Default: \code{XY=NULL}. }
  \item{n.class}{ Number of classes. If \code{n.class=0}, the Sturge equation will be used unless break points are provided. }
  \item{break.pts}{ Vector containing the break points of the distance distribution. Default: \code{break.pts=NULL}. }
  \item{cutoff}{ For the second half of the distance classes, \code{cutoff = TRUE} limits the correlogram to the distance classes that include all points. If \code{cutoff = FALSE}, the correlogram includes all distance classes. }
  \item{r.type}{ Type of correlation in calculation of the Mantel statistic. Default: \code{r.type="pearson"}.  Other choices are \code{r.type="spearman"} and \code{r.type="kendall"}, as in functions \code{\link{cor}} and \code{\link{mantel}}. }
  \item{nperm}{ Number of permutations for the tests of significance. Default: \code{nperm=999}. For large data files, permutation tests are rather slow. } 
  \item{x}{ Output of \code{mantel.correlog}. }
  \item{...}{ Parameters passed to other functions. }
}

\details{
The Mantel tests are done by \code{vegan}'s function \code{\link{mantel}}. Function \code{\link[ape]{mantel.test}} of the \code{ape} package could also have been used. 

The \code{print.mantel.correlog} function prints out the correlogram. See example.
}

\value{ 

  \item{mantel.res }{A table with the distance classes as rows and the class indices, number of distances per class, Mantel statistics, and p-values as columns. }
  \item{n.class }{The number of distance classes. }  
  \item{break.pts }{The break points provided by the user or computed by the program. }  
  \item{call }{The function call. }     
}

\author{ Pierre Legendre, Universite de Montreal }

\references{ 
Legendre, P. & L. Legendre. 1998. Numerical ecology, 2nd English edition. Elsevier Science BV, Amsterdam.

Oden, N. L. & R. R. Sokal. 1986. Directional autocorrelation: an extension of spatial correlograms to two dimensions. Syst. Zool. 35: 608-617.

Sokal, R. R. 1986. Spatial data analysis and historical processes. 29-43 in: E. Diday et al. [eds.] Data analysis and informatics, IV. North-Holland, Amsterdam.
}

\examples{   
# Mite data from "vegan"
data(mite)        
data(mite.xy)  
mite.hel = decostand(mite, "hellinger")
mite.hel.D = dist(mite.hel)

mite.correlog = mantel.correlog(mite.hel.D, XY=mite.xy, nperm=99)
summary(mite.correlog)
mite.correlog
plot(mite.correlog)

mite.correlog2 = mantel.correlog(mite.hel.D, XY=mite.xy, cutoff=FALSE, 
r.type="spearman", nperm=99)
summary(mite.correlog2)
mite.correlog2
plot(mite.correlog2)

}

\keyword{ multivariate }
