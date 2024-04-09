\name{eigenvals}
\Rdversion{1.1}
\alias{eigenvals}
\alias{eigenvals.default}
\alias{eigenvals.prcomp}
\alias{eigenvals.princomp}
\alias{eigenvals.cca}
\alias{eigenvals.wcmdscale}
\alias{eigenvals.pcnm}
\alias{eigenvals.dudi}
\alias{eigenvals.pca}
\alias{eigenvals.pco}
\alias{eigenvals.decorana}
\alias{summary.eigenvals}

\title{
  Extract Eigenvalues from an Ordination Object
}
\description{
  Function extracts eigenvalues from an object that has them. Many
  multivariate methods return such objects. 
}
\usage{
eigenvals(x, ...)
\method{eigenvals}{cca}(x, model = c("all", "unconstrained", "constrained"),
          constrained = NULL, ...)
\method{eigenvals}{decorana}(x, kind = c("additive", "axiswise", "decorana"),
           ...)
\method{summary}{eigenvals}(object, ...)
}

\arguments{
  \item{x}{
    An object from which to extract eigenvalues.
  }
  \item{object}{
    An \code{eigenvals} result object.
  }
  \item{model}{
    Which eigenvalues to return for objects that inherit from class
    \code{"cca"} only.
  }
  \item{constrained}{
    Return only constrained eigenvalues. Deprecated as of vegan
    2.5-0. Use \code{model} instead.
  }
  \item{kind}{
    Kind of eigenvalues returned for \code{\link{decorana}}. Only
    \code{"additive"} eigenvalues can be used for reporting importances
    of components in \code{summary}. \code{"axiswise"} gives the
    non-additive eigenvalues, and \code{"decorana"} the decorana values
    (see \code{\link{decorana}} for details).
  }
  \item{\dots}{
    Other arguments to the functions (usually ignored)
  }
}

\details{
  This is a generic function that has methods for \code{\link{cca}},
  \code{\link{wcmdscale}}, \code{\link{pcnm}}, \code{\link{prcomp}},
  \code{\link{princomp}}, \code{dudi} (of \pkg{ade4}), and 
  \code{pca} and \code{pco} (of 
  \pkg{labdsv}) result objects. The default method also
  extracts eigenvalues if the result looks like being from
  \code{\link{eigen}} or \code{\link{svd}}.  Functions
  \code{\link{prcomp}} and \code{\link{princomp}} contain square roots
  of eigenvalues that all called standard deviations, but
  \code{eigenvals} function returns their squares.  Function
  \code{\link{svd}} contains singular values, but function
  \code{eigenvals} returns their squares. For constrained ordination
  methods \code{\link{cca}}, \code{\link{rda}} and
  \code{\link{capscale}} the function returns the both constrained and
  unconstrained eigenvalues concatenated in one vector, but the partial
  component will be ignored. However, with argument 
  \code{constrained = TRUE} only constrained eigenvalues are returned. 

  The \code{summary} of \code{eigenvals} result returns eigenvalues,
  proportion explained and cumulative proportion explained. The result
  object can have some negative eigenvalues (\code{\link{wcmdscale}},
  \code{\link{dbrda}}, \code{\link{pcnm}}) which correspond to
  imaginary axes of Euclidean mapping of non-Euclidean distances
  (Gower 1985). In these case real axes (corresponding to positive
  eigenvalues) will "explain" proportion >1 of total variation, and
  negative eigenvalues bring the cumulative proportion to
  1. \code{\link{capscale}} will only find the positive eigenvalues
  and only these are used in finding proportions. For
  \code{\link{decorana}} the importances and cumulative proportions
  are only reported for \code{kind = "additive"}, because other
  alternatives do not add up to total inertia of the input data.

}

\value{
  An object of class \code{"eigenvals"}, which is a vector of
  eigenvalues.

  The \code{summary} method returns an object of class
  \code{"summary.eigenvals"}, which is a matrix.
}

\author{
 Jari Oksanen.
}

\references{
  Gower, J. C. (1985). Properties of Euclidean and non-Euclidean
  distance matrices. \emph{Linear Algebra and its Applications} 67,
  81--97.
}

\seealso{
 \code{\link{eigen}}, \code{\link{svd}}, \code{\link{prcomp}},
 \code{\link{princomp}}, \code{\link{cca}}, \code{\link{rda}},
 \code{\link{capscale}}, \code{\link{wcmdscale}},
 \code{\link{cca.object}}. 
}
\examples{
data(varespec)
data(varechem)
mod <- cca(varespec ~ Al + P + K, varechem)
ev <- eigenvals(mod)
ev
summary(ev)

## choose which eignevalues to return
eigenvals(mod, model = "unconstrained")
}
\keyword{ multivariate }

