\encoding{UTF-8}
\name{vegan-deprecated}

\alias{vegan-deprecated}
\alias{as.mlm}
\alias{as.mlm.cca}
\alias{as.mlm.rda}
%------ NOTE:  ../R/vegan-deprecated.R   must be synchronized with this!
\title{Deprecated Functions in vegan package}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)

\description{
  These functions are provided for compatibility with older versions of
  \pkg{vegan} only, and may be defunct as soon as the next release.
}

\usage{
as.mlm(x)
}
\arguments{
  \item{x}{Constrained ordination result. }
}
\details{

  Function \code{as.mlm} refits results of constrained ordination
  (\code{\link{cca}}, \code{\link{rda}}, \code{\link{capscale}}) as a
  multiple response linear model (\code{\link{lm}}). This allows
  finding influence statistics (\code{\link{influence.measures}}).
  It is also possible to derive several other statistics, but most of
  these are biased and misleading, since refitting ignores a major
  component of variation in constrained ordination. This function is
  replaced with a set functions that can find the same statistics
  directly from the ordination result object: see
  \code{\link{hatvalues.cca}}, \code{\link{rstandard.cca}},
  \code{\link{rstudent.cca}}, \code{\link{cooks.distance.cca}},
  \code{\link{vcov.cca}}.

}

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}
