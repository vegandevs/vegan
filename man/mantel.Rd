\name{mantel}
\alias{mantel}
\alias{mantel.partial}
\alias{print.mantel}

\title{Mantel and Partial Mantel Tests for Dissimilarity Matrices }
\description{
  Function \code{mantel}  finds the Mantel statistic as a matrix
  correlation between two dissimilarity matrices, and function
  \code{mantel.partial} finds the partial Mantel statistic as the
  partial matrix correlation between three dissimilarity matricies.  The
  significance of the statistic is evaluated by permuting rows and
  columns of the first dissimilarity matrix. 

}
\usage{
mantel(xdis, ydis, method="pearson", permutations=999, strata)
mantel.partial(xdis, ydis, zdis, method = "pearson", permutations = 999, 
    strata)
}

\arguments{
  \item{xdis, ydis, zdis}{ Dissimilarity matrices or a \code{dist} objects. }
  \item{method}{ Correlation method, as accepted by \code{\link{cor}}:
    \code{"pearson"}, \code{"spearman"} or \code{"kendall"}. }
  \item{permutations}{Number of permutations in assessing significance. }
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
}
\details{
  Mantel statistic is simply a correlation between entries of two
  dissimilarity matrices (some use cross products, but these are linearly
  related).  However, the significance cannot be directly assessed,
  because there are \eqn{N(N-1)/2} entries for just \eqn{N} observations.
  Mantel developed asymptotic test, but here we use permutations of
  \eqn{N} rows and columns of dissimilarity matrix.

  Partial Mantel statistic uses partial correlation
  conditioned on the third matrix. Only the first matrix is permuted so
  that the correlation structure between second and first matrices is
  kept constant. Although \code{mantel.partial} silently accepts other
  methods than \code{"pearson"}, partial correlations will probably be
  wrong with other methods.

  The function uses \code{\link{cor}}, which should accept
  alternatives \code{pearson} for product moment correlations and
  \code{spearman} or \code{kendall} for rank correlations.
}
\value{
  The function returns a list of class \code{mantel} with following
  components: 
  \item{Call }{Function call.}
  \item{method }{Correlation method used, as returned by
    \code{\link[ctest]{cor.test}}.}
  \item{statistic}{The Mantel statistic.}
  \item{signif}{Empirical significance level from permutations.}
  \item{perm}{A vector of permuted values.}
  \item{permutations}{Number of permutations.}
}
\references{ The test is due to Mantel, of course, but the
  current implementation is based on Legendre and Legendre.

  Legendre, P. and Legendre, L. (1998) \emph{Numerical Ecology}. 2nd English
  Edition. Elsevier.
  
}

\note{
  Legendre & Legendre (1998) say that partial Mantel correlations 
  often are difficult to interpet. 
  }

\author{Jari Oksanen }


\seealso{\code{\link{cor}} for correlation coefficients,
  \code{\link{protest}} (``Procrustes test'') for an alternative with
  ordination diagrams, \code{\link{anosim}}  and \code{\link{mrpp}}
  for comparing dissimilarities against 
  classification.  For dissimilarity matrices, see \code{\link{vegdist}}
  or \code{\link{dist}}.  See \code{\link{bioenv}} for selecting
  environmental variables. }

\examples{
## Is vegetation related to environment?
data(varespec)
data(varechem)
veg.dist <- vegdist(varespec) # Bray-Curtis
env.dist <- vegdist(scale(varechem), "euclid")
mantel(veg.dist, env.dist)
mantel(veg.dist, env.dist, method="spear")
}
\keyword{ multivariate }
\keyword{ htest }
