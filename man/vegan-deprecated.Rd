\name{vegan-deprecated}

\alias{orditkplot}
\alias{as.mcmc.oecosimu}
\alias{as.mcmc.permat}
\alias{vegan-deprecated}

%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Deprecated Functions in vegan package}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)

\description{
  These functions are provided for compatibility with older versions of
  \pkg{vegan} only, and may be defunct as soon as the next release.
}

\usage{

## moved to vegan3d package: install from CRAN
orditkplot(...)
## use toCoda instead
as.mcmc.oecosimu(x)
as.mcmc.permat(x)
}

\arguments{

  \item{x}{object to be transformed.}
  \item{\dots}{Other arguments passed to the functions.}

}

\details{

  \code{orditkplot} was moved to \CRANpkg{vegan3d} (version
  1.3-0). Install that package from CRAN and use in the old way.

  \code{as.mcmc} functions were replaced with \code{\link{toCoda}}.
}

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}
