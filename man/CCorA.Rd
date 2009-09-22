\name{CCorA}
\alias{CCorA}
\alias{print.CCorA}
\alias{biplot.CCorA}
\concept{ordination}

\title{Canonical Correlation Analysis}

\description{Canonical correlation analysis, following Brian McArdle's
unpublished graduate course notes, plus improvements to allow the
calculations in the case of very sparse and collinear matrices. }

\usage{
CCorA(Y, X, stand.Y=FALSE, stand.X=FALSE, nperm = 0, ...)

\method{biplot}{CCorA}(x, xlabs, which = 1:2, ...)
}

\arguments{
  \item{Y}{ left matrix. }
  \item{X}{ right matrix. }
  \item{stand.Y}{ logical; should \code{Y} be standardized? }
  \item{stand.X}{ logical; should \code{X} be standardized? }
  \item{nperm}{ numeric; Number of permutations to evaluate the
    significance of Pillai's trace}
  \item{x}{\code{CCoaR} result object}
  \item{xlabs}{Row labels. The default is to use row names, \code{NULL}
    uses row numbers instead, and \code{NA} suppresses plotting row names
    completely}
  \item{which}{ \code{1} plots \code{Y} results, and
    \code{2} plots \code{X1} results }
  \item{\dots}{Other arguments passed to functions. \code{biplot.CCorA}
    passes graphical arguments to \code{\link{biplot}} and
    \code{\link{biplot.default}}, \code{CCorA} currently ignores extra 
    arguments.} 
}

\details{
  Canonical correlation analysis (Hotelling 1936) seeks linear
  combinations of the variables of \code{Y} that are maximally
  correlated to linear combinations of the variables of \code{X}. The
  analysis estimates the relationships and displays them in graphs.

  Algorithmic notes:
  \enumerate{
    \item
    All data matrices are replaced by their PCA object scores, computed
    by SVD.
    \item
    The blunt approach would be to read the three matrices, compute the
    covariance matrices, then the matrix
    \code{S12 \%*\% inv(S22) \%*\% t(S12) \%*\% inv(S11)}.
    Its trace is Pillai's trace statistic. 
    \item
    This approach may fail, however, when there is heavy multicollinearity
    in very sparse data matrices, as it is the case in 4th-corner inflated
    data matrices for example. The safe approach is to replace all data
    matrices by their PCA object scores.
    \item
    Inversion by \code{\link{solve}} is avoided. Computation of inverses
    is done by \acronym{SVD}  (\code{\link{svd}}) in most cases.
    \item
    Regression by \acronym{OLS} is also avoided. Regression residuals are
    computed by \acronym{QR} decomposition (\code{\link{qr}}).
  }

The \code{biplot} function can produce two biplots, each for the left
matrix and right matrix solutions. The function passes all arguments to
\code{\link{biplot.default}}, and you should consult its help page for
configuring biplots.
}

\value{
  Function \code{CCorA} returns a list containing the following components:
  \item{ Pillai }{ Pillai's trace statistic = sum of canonical
    eigenvalues. } 
  \item{ EigenValues }{ Canonical eigenvalues. They are the squares of the
    canonical correlations. }
  \item{ CanCorr }{ Canonical correlations. }
  \item{ Mat.ranks }{ Ranks of matrices \code{Y} and X1 (possibly after
    controlling for X2). }
  \item{ RDA.Rsquares }{ Bimultivariate redundancy coefficients
    (R-squares) of RDAs of Y|X1 and X1|Y. } 
  \item{ RDA.adj.Rsq }{ \code{RDA.Rsquares} adjusted for n and number of
    explanatory variables. }
  \item{ AA }{ Scores of Y variables in Y biplot. }
  \item{ BB }{ Scores of X1 variables in X1 biplot. }
  \item{ Cy }{ Object scores in Y biplot. }
  \item{ Cx }{ Object scores in X1 biplot. }
}

\references{ 
  Hotelling, H. 1936. Relations between two sets of
  variates. \emph{Biometrika} \strong{28}: 321-377.
}

\author{ Pierre Legendre, Departement de Sciences Biologiques,
  Universite de Montreal. Implemented in \pkg{vegan} with the help of
  Jari Oksanen. }

\examples{
# Example using random numbers
mat1 <- matrix(rnorm(60),20,3)
mat2 <- matrix(rnorm(100),20,5)
CCorA(mat1, mat2)

# Example using intercountry life-cycle savings data, 50 countries
data(LifeCycleSavings)
pop <- LifeCycleSavings[, 2:3]
oec <- LifeCycleSavings[, -(2:3)]
out <- CCorA(pop, oec)
out
biplot(out, xlabs = NA)
}

\keyword{ multivariate }
