\name{rankindex}
\alias{rankindex}

\title{Compares Dissimilarity Indices for Gradient Detection }
\description{
  Rank correlations between dissimilarity indices
  and gradient separation.
}
\usage{
rankindex(grad, veg, indices = c("euc", "man", "gow", "bra", "kul"),
          stepacross = FALSE, method = "spearman", 
	  metric = c("euclidean", "mahalanobis", "manhattan", "gower"),
	  ...)
}

\arguments{
  \item{grad}{The gradient variable or matrix. }
  \item{veg}{The community data matrix. }
  \item{indices}{Dissimilarity indices compared, partial matches to
    alternatives in \code{\link{vegdist}}.
    Alternatively, it can be a (named) list of functions returning
    objects of class 'dist'.}
  \item{stepacross}{Use \code{\link{stepacross}} to find
    a shorter path dissimilarity. The dissimilarities for site pairs
    with no shared species are set \code{NA} using
    \code{\link{no.shared}} so that indices with no fixed
    upper limit can also be analysed.}
  \item{method}{Correlation method used.}
  \item{metric}{Metric to evaluate the gradient separation. See Details.}
  \item{...}{Other parameters to \code{\link{stepacross}}.}
}
\details{
  A good dissimilarity index for multidimensional scaling should have
  a high rank-order similarity with gradient separation.  The function
  compares most indices in \code{\link{vegdist}} against gradient
  separation using rank correlation coefficients in
  \code{\link{cor}}. The gradient separation between each point is
  assessed using given \code{metric}. The default is to use Euclidean
  distance of continuous variables scaled to unit variance, or to use
  Gower metric for mixed data using function
  \code{\link[cluster]{daisy}} when \code{grad} has factors. The other
  alternatives are Mahalanabis distances which are based on
  \code{grad} matrix scaled so that columns are orthogonal
  (uncorrelated) and have unit variance, or Manhattan distances of
  \code{grad} variables scaled to unit range.

  The \code{indices} argument can accept any dissimilarity 
  indices besides the ones calculated by the 
  \code{\link{vegdist}} function. For this, the argument value
  should be a (possibly named) list of functions.
  Each function must return a valid 'dist' object with dissimilarities,
  similarities are not accepted and should be converted into dissimilarities
  beforehand.
}
\value{
  Returns a named vector of rank correlations.
}
\references{ Faith, F.P., Minchin, P.R. and Belbin,
  L. (1987).  Compositional dissimilarity as a robust measure of
    ecological distance. \emph{Vegetatio} 69, 57-68. }
\author{Jari Oksanen, with additions from Peter Solymos}
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

  The functions returning dissimilarity objects should be self contained,
  because the \code{...} argument passes additional parameters
  to \code{\link{stepacross}} and not to the functions supplied
  via the \code{indices} argument.
}

\seealso{\code{\link{vegdist}}, \code{\link{stepacross}},
  \code{\link{no.shared}}, \code{\link{monoMDS}},
    \code{\link{cor}}, \code{\link{Machine}}, and for
    alternatives \code{\link{anosim}}, \code{\link{mantel}} and
    \code{\link{protest}}. }
  
\examples{
data(varespec)
data(varechem)
## The variables are automatically scaled
rankindex(varechem, varespec)
rankindex(varechem, wisconsin(varespec))
## Using non vegdist indices as functions
funs <- list(Manhattan=function(x) dist(x, "manhattan"),
    Gower=function(x) cluster:::daisy(x, "gower"),
    Ochiai=function(x) designdist(x, "1-J/sqrt(A*B)"))
rankindex(scale(varechem), varespec, funs)
}
\keyword{ multivariate }
