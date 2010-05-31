\encoding{UTF-8}
\name{mantel.correlog}
\alias{mantel.correlog}
\alias{print.mantel.correlog}
\alias{plot.mantel.correlog}
\title{ Mantel Correlogram }

\description{
  Function \code{mantel.correlog} computes a multivariate
  Mantel correlogram. Proposed by Sokal (1986) and Oden and Sokal
  (1986), the method is also described in Legendre and Legendre (1998,
  pp. 736-738).
}

\usage{
mantel.correlog(D.eco, D.geo=NULL, XY=NULL, n.class=0, break.pts=NULL, 
cutoff=TRUE, r.type="pearson", nperm=999, mult="holm", progressive=TRUE)
\method{plot}{mantel.correlog}(x, alpha=0.05, ...)
}

\arguments{
  \item{D.eco}{ An ecological distance matrix, with class
  either \code{dist} or \code{matrix}. }
  
  \item{D.geo}{ A geographic distance matrix, with class either
  \code{dist} or \code{matrix}. Provide either \code{D.geo} or
  \code{XY}. Default: \code{D.geo=NULL}. }

  \item{XY}{ A file of Cartesian geographic coordinates of the
  points. Default: \code{XY=NULL}. }

  \item{n.class}{ Number of classes. If \code{n.class=0}, the Sturges
  equation will be used unless break points are provided. }

  \item{break.pts}{ Vector containing the break points of the distance
  distribution. Provide (n.class+1) breakpoints, that is, a list with
  a beginning and an ending point. Default: \code{break.pts=NULL}. }

  \item{cutoff}{ For the second half of the distance classes,
  \code{cutoff = TRUE} limits the correlogram to the distance classes
  that include all points. If \code{cutoff = FALSE}, the correlogram
  includes all distance classes. }

  \item{r.type}{ Type of correlation in calculation of the Mantel
  statistic. Default: \code{r.type="pearson"}.  Other choices are
  \code{r.type="spearman"} and \code{r.type="kendall"}, as in functions
  \code{\link{cor}} and \code{\link{mantel}}. }

  \item{nperm}{ Number of permutations for the tests of
  significance. Default: \code{nperm=999}. For large data files,
  permutation tests are rather slow. }

  \item{mult}{ Correct P-values for multiple testing. The correction
  methods are \code{"holm"} (default), \code{"hochberg"},
  \code{"sidak"}, and other methods available in the
  \code{\link{p.adjust}} function: \code{"bonferroni"} (best known, but
  not recommended because it is overly conservative), \code{"hommel"},
  \code{"BH"}, \code{"BY"}, \code{"fdr"}, and \code{"none"}. }

  \item{progressive}{ Default: \code{progressive=TRUE} for progressive
  correction of multiple-testing, as described in Legendre and Legendre
  (1998, p. 721). Test of the first distance class: no correction;
  second distance class: correct for 2 simultaneous tests; distance
  class k: correct for k simultaneous tests. \code{progressive=FALSE}:
  correct all tests for \code{n.class} simultaneous tests. }

  \item{x}{ Output of \code{mantel.correlog}. }

  \item{alpha}{ Significance level for the points drawn with black
  symbols in the correlogram. Default: \code{alpha=0.05}. }

  \item{...}{ Other parameters passed from other functions. }
}

\details{ A correlogram is a graph in which spatial correlation values
  are plotted, on the ordinate, as a function of the geographic distance
  classes among the study sites along the abscissa. In a Mantel
  correlogram, a Mantel correlation (Mantel 1967) is computed between a
  multivariate (e.g. multi-species) distance matrix of the user's choice
  and a design matrix representing each of the geographic distance
  classes in turn. The Mantel statistic is tested through a
  permutational Mantel test performed by \code{vegan}'s
  \code{\link{mantel}} function.

  When a correction for multiple testing is applied, more permutations
  are necessary than in the no-correction case, to obtain significant
  p-values in the higher correlogram classes.

  The \code{print.mantel.correlog} function prints out the
  correlogram. See examples.  }

\value{ 

  \item{mantel.res }{A table with the distance classes as rows and the
  class indices, number of distances per class, Mantel statistics
  (computed using Pearson's r, Spearman's r, or Kendall's tau), and
  p-values as columns. A positive Mantel statistic indicates positive
  spatial correlation. An additional column with p-values corrected for
  multiple testing is added unless \code{mult="none"}. }

  \item{n.class }{The n umber of distance classes. }
  
  \item{break.pts }{The break points provided by the user or computed by
    the program. }

  \item{mult }{The name of the correction for multiple testing. No
    correction: \code{mult="none"}. }  

  \item{progressive }{A logical (\code{TRUE}, \code{FALSE}) value
  indicating whether or not a progressive correction for multiple
  testing was requested. } 

  \item{n.tests }{The number of distance classes for which Mantel
  tests have been computed and tested for significance. }

  \item{call }{The function call. }  
}

\author{ Pierre Legendre, Université de Montréal }

\references{

  Legendre, P. and L. Legendre. 1998. Numerical ecology, 2nd English
  edition. Elsevier Science BV, Amsterdam.

  Mantel, N. 1967. The detection of disease clustering and a generalized
  regression approach. Cancer Res. 27: 209-220.

  Oden, N. L. and R. R. Sokal. 1986. Directional autocorrelation: an
  extension of spatial correlograms to two dimensions. Syst. Zool. 35:
  608-617.

  Sokal, R. R. 1986. Spatial data analysis and historical
  processes. 29-43 in: E. Diday et al. [eds.] Data analysis and
  informatics, IV. North-Holland, Amsterdam.
  
  Sturges, H. A. 1926. The choice of a class interval. Journal of the 
  American Statistical Association 21: 65–66.  }

\examples{   
# Mite data available in "vegan"
data(mite)        
data(mite.xy)  
mite.hel <- decostand(mite, "hellinger")

# Detrend the species data by regression on the site coordinates
mite.hel.resid <- resid(lm(as.matrix(mite.hel) ~ ., data=mite.xy))

# Compute the detrended species distance matrix
mite.hel.D = dist(mite.hel.resid)

# Compute Mantel correlogram with cutoff, Pearson statistic
mite.correlog = mantel.correlog(mite.hel.D, XY=mite.xy, nperm=49)
summary(mite.correlog)
mite.correlog   
# or: print(mite.correlog)
# or: print.mantel.correlog(mite.correlog)
plot(mite.correlog)

# Compute Mantel correlogram without cutoff, Spearman statistic
mite.correlog2 = mantel.correlog(mite.hel.D, XY=mite.xy, cutoff=FALSE, 
r.type="spearman", nperm=49)
summary(mite.correlog2)
mite.correlog2
plot(mite.correlog2)

# NOTE: 'nperm' argument usually needs to be larger than 49.
# It was set to this low value for demonstration purposes.

}

\keyword{ multivariate }