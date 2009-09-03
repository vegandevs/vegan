\name{rankindex}
\alias{rankindex}

\title{Compares Dissimilarity Indices for Gradient Detection }
\description{
  Rank correlations between dissimilarity indices
  and gradient separation.
}
\usage{
rankindex(grad, veg, indices = c("euc", "man", "gow", "bra", "kul"),
          stepacross = FALSE, method = "spearman", ...)
}

\arguments{
  \item{grad}{The gradient variable or matrix. }
  \item{veg}{The community data matrix. }
  \item{indices}{Dissimilarity indices compared, partial matches to
    alternatives in \code{\link{vegdist}}.}
  \item{stepacross}{Use \code{\link{stepacross}} to find
    a shorter path dissimilarity. The dissimilarities for site pairs
    with no shared species are set \code{NA} using
    \code{\link{no.shared}} so that indices with no fixed
    upper limit can also be analysed.}
  \item{method}{Correlation method used. }
  \item{...}{Other parameters to \code{\link{stepacross}}.}
}
\details{
  A good dissimilarity index for multidimensional scaling 
  should have a high rank-order similarity with gradient separation.
  The function compares most indices in \code{\link{vegdist}} against
  gradient separation using rank correlation coefficients in
  \code{\link{cor.test}}. The gradient separation between each
  point is assessed as Euclidean distance for continuous variables, and
  as Gower metric for mixed data using function
  \code{\link[cluster]{daisy}} when \code{grad} has factors.
}
\value{
  Returns a named vector of rank correlations.
}
\references{ Faith, F.P., Minchin, P.R. and Belbin,
  L. (1987).  Compositional dissimilarity as a robust measure of
    ecological distance. \emph{Vegetatio} 69, 57-68. }
\author{Jari Oksanen }
\note{
  There are several problems in using rank correlation coefficients.
  Typically there are very many ties when \eqn{n(n-1)/2} gradient
  separation values are derived from just \eqn{n} observations.
  Due to floating point arithmetics, many tied values differ by
  machine epsilon and are arbitrarily ranked differently by
  \code{\link{rank}} used in \code{\link{cor.test}}.  Two indices
  which are identical with certain
  transformation or standardization may differ slightly
  (magnitude \eqn{10^{-15}}) and this may lead into third or fourth decimal
  instability in rank correlations.  Small differences in rank
  correlations should not be taken too seriously.  Probably this method
  should be replaced with a sounder method, but I do not yet know
  which\ldots  You may experiment with \code{\link{mantel}},
  \code{\link{anosim}} or even \code{\link{protest}}.

  Earlier version of this function used \code{method = "kendall"}, but
  that is far too slow in large data sets.
}

\seealso{\code{\link{vegdist}}, \code{\link{stepacross}},
  \code{\link{no.shared}}, \code{\link[MASS]{isoMDS}},
    \code{\link{cor}}, \code{\link{Machine}}, and for
    alternatives \code{\link{anosim}}, \code{\link{mantel}} and
    \code{\link{protest}}. }
  
\examples{
data(varespec)
data(varechem)
## The next scales all environmental variables to unit variance.
## Some would use PCA transformation.
rankindex(scale(varechem), varespec)
rankindex(scale(varechem), wisconsin(varespec))
}
\keyword{ multivariate }
