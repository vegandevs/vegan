\name{vegan-deprecated}
\alias{metaMDSrotate}

\alias{vegan-deprecated}
%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Deprecated Functions in vegan package}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  These functions are provided for compatibility with older versions of
  \pkg{vegan} only, and may be defunct as soon as the next release.
}
\usage{
metaMDSrotate(object, vec, na.rm = FALSE, ...)
}

\arguments{
 \item{object}{A result object from \code{metaMDS}.}
 \item{vec}{A continuous site variable (vector).}
 \item{na.rm}{Remove missing values from continuous variable \code{vec}.}
  \item{\dots}{Other parameters passed to functions.}
}

\details{
  %-- Either:
  %There are currently no deprecated functions in this package.
  %-- Or:
  %% explain *why* it's deprecated, and \code{\link{..}} to new
  
  Function \code{metaMDSrotate} is replaced with
  \code{\link{MDSrotate}} which can handle \code{\link{monoMDS}}
  results in addition to \code{\link{metaMDS}}.

}
\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}
