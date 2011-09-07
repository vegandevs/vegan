\name{vegan-defunct}
%--- The following functions were moved to the 'permute' package and
%    removed from vegan, but here we document only those that were
%    renamed and are not documented in 'permute'
 
%\alias{permCheck}
%\alias{numPerms}
%\alias{print.permCheck}
%\alias{print.summary.permCheck}
%\alias{summary.permCheck}
\alias{getNumObs}
%\alias{getNumObs.default}
%\alias{getNumObs.integer}
%\alias{getNumObs.numeric}
%\alias{allPerms}
%\alias{print.allPerms}
%\alias{summary.allPerms}
%\alias{print.summary.allPerms}
%\alias{permuplot}
\alias{permuted.index2}
%\alias{permControl}
%\alias{print.permControl}
%\alias{permute}

\alias{vegan-defunct}
%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Defunct Functions in Package \pkg{vegan}}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  The functions or variables listed here are no longer part of \pkg{vegan} as
  they are no longer needed.
}
\usage{
%-- Removed from vegan 2.0-0: now in package permute
getNumObs(object, \dots)
\method{getNumObs}{default}(object, \dots)
\method{getNumObs}{numeric}(object, \dots)
\method{getNumObs}{integer}(object, \dots)

permuted.index2(n, control = permControl())
}

\details{ The \dQuote{new} permutation functions were moved to the
  \pkg{permute} package, and they are documented there.  The
  \pkg{permute} package replaces \code{permuted.index2} with
  \code{\link[permute]{shuffle}} and \code{getNumObs} with its specific
  \code{\link[permute]{nobs-methods}}.
}

\seealso{
  \code{\link{Defunct}}, \code{\link{vegan-deprecated}}
}

\keyword{internal}
