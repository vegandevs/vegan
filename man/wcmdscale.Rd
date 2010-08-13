% Based on:
% File src/library/stats/man/cmdscale.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{wcmdscale}
\alias{wcmdscale}
\title{Weighted Classical (Metric) Multidimensional Scaling}
\usage{
wcmdscale(d, k, eig = FALSE, add = FALSE, x.ret = FALSE, w)
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
  \item{add}{logical indicating if an additive constant \eqn{c*} should
    be computed, and added to the non-diagonal dissimilarities such that
    all \eqn{n-1} eigenvalues are non-negative. \strong{Not implemented}. }
  \item{x.ret}{indicates whether the doubly centred symmetric distance
    matrix should be returned.}
  \item{w}{Weights of points.}
}
\details{
  Function \code{wcmdscale} is based on function
  \code{\link{cmdscale}} (package \pkg{stats} of base \R), but it uses
  point weights. Points with high weights will have a stronger
  influence on the result than those with low weights. Setting equal
  weights \code{w = 1} will give ordinary multidimensional scaling.
}
\value{
  If \code{eig = FALSE} and \code{x.ret = FALSE} (default), a matrix
  with \code{k} columns whose rows give the coordinates of the points
  chosen to represent the dissimilarities.

  Otherwise, an object of class \code{wcmdscale} list containing the
  following components. 
  \item{points}{a matrix with \code{k} columns whose rows give the
    coordinates of the points chosen to represent the dissimilarities.}
  \item{eig}{the \eqn{n-1} eigenvalues computed during the scaling process if
    \code{eig} is true.}
  \item{x}{the doubly centred and weighted distance matrix if \code{x.ret} is true.}
  \item{weights}{Weights.}
  \item{negaxes}{A matrix of scores for axes with negative eigenvalues scaled 
     by the absolute eigenvalues similarly as \code{points}. This is \code{NULL}
     if there are no negative eigenvalues or \code{k} was specified, and would not     include negative eigenvalues.}
}
\references{
  Gower, J. C. (1966)  
  Some distance properties of latent root and vector 
  methods used in multivariate analysis.  
  \emph{Biometrika} \bold{53}, 325--328.

  Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979).  Chapter 14 of
  \emph{Multivariate Analysis}, London: Academic Press.
}
\seealso{
  \code{\link{cmdscale}}.
  Also \code{\link[MASS]{isoMDS}} and \code{\link[MASS]{sammon}}
  in package \pkg{MASS}.
}
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
procrustes(ord, ca, choices=1:19, scaling = 1)
plot(procrustes(ord, ca, choices=1:2, scaling=1))
## Reconstruction of non-Euclidean distances with negative eigenvalues
d <- vegdist(dune)
ord <- wcmdscale(d, eig = TRUE)
## Only positive eigenvalues:
cor(d, dist(ord$points))
## Correction with negative eigenvalues:
cor(d, sqrt(dist(ord$points)^2 - dist(ord$negaxes)^2))
}
\keyword{multivariate}
