\name{vegan-deprecated}

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

## use toCoda instead
as.mcmc.oecosimu(x)
as.mcmc.permat(x)
}

\arguments{

  \item{x}{object to be transformed.}

}

\details{

  \code{as.mcmc} functions were replaced with \code{\link{toCoda}}.
}

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}
