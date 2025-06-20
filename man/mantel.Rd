\name{mantel}
\alias{mantel}
\alias{mantel.partial}
\alias{summary.mantel}

\title{Mantel and Partial Mantel Tests for Dissimilarity Matrices }
\description{
  Function \code{mantel}  finds the Mantel statistic as a matrix
  correlation between two dissimilarity matrices, and function
  \code{mantel.partial} finds the partial Mantel statistic as the
  partial matrix correlation between three dissimilarity matrices.  The
  significance of the statistic is evaluated by permuting rows and
  columns of the first dissimilarity matrix. Test is one-sided and only
  tests that distances are positively correlated.

}
\usage{
mantel(xdis, ydis, method="pearson", permutations=999, strata = NULL,
    na.rm = FALSE, parallel = getOption("mc.cores"))
mantel.partial(xdis, ydis, zdis, method = "pearson", permutations = 999, 
    strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))
\method{summary}{mantel}(object, ...)
}

\arguments{

  \item{xdis, ydis, zdis}{Distance object of class \code{"dist"} or
    symmetric square matrices of distances. Only the lower triangle of
    square matrices is used.  The first object \code{xdis} will be
    permuted in permutation tests. }
  \item{method}{ Correlation method, as accepted by \code{\link{cor}}:
    \code{"pearson"}, \code{"spearman"} or \code{"kendall"}. }
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
  \item{na.rm}{Remove missing values in calculation of Mantel
    correlation. Use this option with care: Permutation tests can
    be biased, in particular if two matrices had missing values in
    matching positions.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}
  \item{object}{Result object.}
  \item{\dots}{Arguments passed to \code{\link{summary.permustats}}
    These include \code{alternative} to select the sidedness of the
    test.}

}

\details{
  Mantel statistic is simply a correlation between entries of two
  dissimilarity matrices (some use cross products, but these are
  linearly related).  However, the significance cannot be directly
  assessed, because there are \eqn{N(N-1)/2} entries for just \eqn{N}
  observations.  Mantel developed asymptotic test, but here we use
  permutations of \eqn{N} rows and columns of dissimilarity
  matrix. Only the first matrix (\code{xdist}) will be permuted, and
  the second is kept constant. See \code{\link{permutations}} for
  additional details on permutation tests in Vegan.

  Partial Mantel statistic uses partial correlation
  conditioned on the third matrix. Only the first matrix is permuted so
  that the correlation structure between second and first matrices is
  kept constant. Although \code{mantel.partial} silently accepts other
  methods than \code{"pearson"}, partial correlations will probably be
  wrong with other methods.

  Borcard & Legendre (2012) warn against using partial Mantel test and
  recommend instead Mantel correlogram
  (\code{\link{mantel.correlog}}).

  The function uses \code{\link{cor}}, which should accept
  alternatives \code{pearson} for product moment correlations and
  \code{spearman} or \code{kendall} for rank correlations.
}
\value{
  The function returns a list of class \code{mantel} with following
  components: 
  \item{Call }{Function call.}
  \item{method }{Correlation method used, as returned by
    \code{\link{cor.test}}.}
  \item{statistic}{The Mantel statistic.}
  \item{signif}{Empirical significance level from permutations.}
  \item{perm}{A vector of permuted values. The distribution of
    permuted values can be inspected with \code{\link{permustats}} 
    function.}
  \item{permutations}{Number of permutations.}
  \item{control}{A list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}.}
}
\references{

  Borcard, D. & Legendre, P. (2012) Is the Mantel correlogram powerful
  enough to be useful in ecological analysis? A simulation
  study. \emph{Ecology} 93: 1473-1481.

  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd English
  Edition. Elsevier.
  
}


\author{Jari Oksanen }


\seealso{ \code{\link{mantel.correlog}}. }

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
