\name{vegan-deprecated}

\alias{adonis}
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
## use adonis2 instead
adonis(formula, data, permutations = 999, method = "bray",
    strata = NULL, contr.unordered = "contr.sum",
    contr.ordered = "contr.poly", parallel = getOption("mc.cores"), ...)
## use toCoda instead
as.mcmc.oecosimu(x)
as.mcmc.permat(x)
}

\arguments{
  \item{formula, data, permutations, method, parallel, \dots}{See
    \code{\link{adonis2}}.}
  \item{strata}{groups (strata) within which to constrain permutations.}
  \item{contr.unordered, contr.ordered}{contrasts used for design matrix.}
  \item{x}{object to be tranformed.}
}

\details{

  \code{\link{adonis2}} replaces \code{adonis} with extended
  functionality and completely new internal design. The shared arguments
  of \code{adonis} are similar as in \code{\link{adonis2}}, but
  arguments \code{contr.unordered} and \code{contr.ordered} can set the
  contrasts within \code{adonis}.

}

\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}
