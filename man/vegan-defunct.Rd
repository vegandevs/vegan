\name{vegan-defunct}
%--- The following functions were moved to the 'permute' package and
%    removed from vegan, but here we document only those that were
%    renamed and are not documented in 'permute'
 
\alias{metaMDSrotate}

\alias{vegan-defunct}
%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Defunct Functions in Package \pkg{vegan}}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  The functions or variables listed here are no longer part of \pkg{vegan} as
  they are no longer needed.
}
\usage{
metaMDSrotate(object, vec, na.rm = FALSE, ...)
}

\details{ 
  Function \code{metaMDSrotate} is replaced with
  \code{\link{MDSrotate}} which can handle \code{\link{monoMDS}}
  results in addition to \code{\link{metaMDS}}.
}

\seealso{
  \code{\link{Defunct}}, \code{\link{vegan-deprecated}}
}

\keyword{internal}
