% Based on:
% File src/library/stats/man/cmdscale.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{wcmdscale}
\alias{wcmdscale}
\alias{scores.wcmdscale}
\alias{plot.wcmdscale}

\title{Weighted Classical (Metric) Multidimensional Scaling}
\usage{
wcmdscale(d, k, eig = FALSE, add = FALSE, x.ret = FALSE, w)
\method{plot}{wcmdscale}(x, choices = c(1, 2), type = "t", ...)
\method{scores}{wcmdscale}(x, choices = NA, ...)
}
\description{
  Weighted classical multidimensional scaling,
  also known as weighted \emph{principal coordinates analysis}.
}
\arguments{
  \item{d}{a distance structure such as that returned by \code{dist}
    or a full symmetric matrix containing the dissimilarities.}
  \item{k}{the dimension of the space which the data are to be
    represented in; must be in \eqn{\{1,2,\ldots,n-1\}}{{1,2,\ldots,n-1}}.
    If missing, all dimensions with above zero eigenvalue.}
  \item{eig}{indicates whether eigenvalues should be returned.}
  \item{add}{an additive constant \eqn{c} is added to the non-diagonal
    dissimilarities such that all \eqn{n-1} eigenvalues are
    non-negative. Alternatives are \code{"lingoes"} (default, also
    used with \code{TRUE}) and \code{"cailliez"} (which is the only
    alternative in \code{\link{cmdscale}}). See Legendre & Anderson
    (1999).}
  \item{x.ret}{indicates whether the doubly centred symmetric distance
    matrix should be returned.}
  \item{w}{Weights of points.}
  \item{x}{The \code{wcmdscale} result object when the function was
    called with options \code{eig = TRUE} or \code{x.ret = TRUE} (See
    Details).}
  \item{choices}{Axes to be returned; \code{NA} returns all real axes.}
  \item{type}{Type of graph which may be \code{"t"}ext, \code{"p"}oints
    or \code{"n"}one.}
  \item{\dots}{Other arguments passed to graphical functions.}
}
\details{
  Function \code{wcmdscale} is based on function
  \code{\link{cmdscale}} (package \pkg{stats} of base \R), but it uses
  point weights. Points with high weights will have a stronger
  influence on the result than those with low weights. Setting equal
  weights \code{w = 1} will give ordinary multidimensional scaling.

  With default options, the function returns only a matrix of scores
  scaled by eigenvalues for all real axes. If the function is called
  with \code{eig = TRUE} or \code{x.ret = TRUE}, the function returns
  an object of class \code{"wcmdscale"} with \code{print},
  \code{plot}, \code{scores}, \code{\link{eigenvals}} and
  \code{\link{stressplot}} methods, and described in section Value.

  The method is Euclidean, and with non-Euclidean dissimilarities some
  eigenvalues can be negative. If this disturbs you, this can be
  avoided by adding a constant to non-diagonal dissimilarities making
  all eigenvalues non-negative. The function implements methods
  discussed by Legendre & Anderson (1999): The method of Lingoes
  (\code{add="lingoes"}) adds the constant \eqn{c} to squared
  dissimilarities \eqn{d} using \eqn{\sqrt{d^2 + 2 c}}{sqrt(d^2 + 2*c)}
  and the method of Cailliez (\code{add="cailliez"}) to
  dissimilarities using \eqn{d + c}. Legendre & Anderson (1999)
  recommend the method of Lingoes, and base \R{} function
  \code{\link{cmdscale}} implements the method of Cailliez.
}

\value{ If \code{eig = FALSE} and \code{x.ret = FALSE} (default), a
  matrix with \code{k} columns whose rows give the coordinates of
  points corresponding to positive eigenvalues.  Otherwise, an object
  of class \code{wcmdscale} containing the components that are mostly
  similar as in \code{\link{cmdscale}}:

  \item{points}{a matrix with \code{k} columns whose rows give the
    coordinates of the points chosen to represent the
    dissimilarities.}
  \item{eig}{the \eqn{n-1} eigenvalues computed during the scaling
    process if \code{eig} is true.}
  \item{x}{the doubly centred and weighted distance matrix if
    \code{x.ret} is true.}
  \item{ac, add}{additive constant and adjustment method used to avoid
    negative eigenvalues. These are \code{NA} and \code{FALSE} if no
    adjustment was done.}
  \item{GOF}{Goodness of fit statistics for \code{k} axes. The first
    value is based on the sum of absolute values of all eigenvalues,
    and the second value is based on the sum of positive eigenvalues}
  \item{weights}{Weights.}
  \item{negaxes}{A matrix of scores for axes with negative eigenvalues
     scaled by the absolute eigenvalues similarly as
     \code{points}. This is \code{NULL} if there are no negative
     eigenvalues or \code{k} was specified, and would not include
     negative eigenvalues.}
  \item{call}{Function call.}
}

\references{
  Gower, J. C. (1966)
  Some distance properties of latent root and vector
  methods used in multivariate analysis.
  \emph{Biometrika} \bold{53}, 325--328.

  Legendre, P. & Anderson, M. J. (1999). Distance-based redundancy
    analysis: testing multispecies responses in multifactorial
    ecological experiments. \emph{Ecology} \bold{69}, 1--24.

  Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979).  Chapter 14 of
  \emph{Multivariate Analysis}, London: Academic Press.
}

\seealso{The function is modelled after \code{\link{cmdscale}}, but adds
  weights (hence name) and handles negative eigenvalues differently.
  \code{\link{eigenvals.wcmdscale}} and
  \code{\link{stressplot.wcmdscale}} are some specific methods. Other
  multidimensional scaling methods are \code{\link{monoMDS}}, and
  \code{\link[MASS]{isoMDS}} and \code{\link[MASS]{sammon}} in package
  \pkg{MASS}.  }

\examples{
## Correspondence analysis as a weighted principal coordinates
## analysis of Euclidean distances of Chi-square transformed data
data(dune)
rs <- rowSums(dune)/sum(dune)
d <- dist(decostand(dune, "chi"))
ord <- wcmdscale(d, w = rs, eig = TRUE)
## Ordinary CA
ca <- cca(dune)
## Eigevalues are numerically similar
ca$CA$eig - ord$eig
## Configurations are similar when site scores are scaled by
## eigenvalues in CA
procrustes(ord, ca, choices=1:19, scaling = "sites")
plot(procrustes(ord, ca, choices=1:2, scaling="sites"))
## Reconstruction of non-Euclidean distances with negative eigenvalues
d <- vegdist(dune)
ord <- wcmdscale(d, eig = TRUE)
## Only positive eigenvalues:
cor(d, dist(ord$points))
## Correction with negative eigenvalues:
cor(d, sqrt(dist(ord$points)^2 - dist(ord$negaxes)^2))
}
\keyword{multivariate}
