\encoding{UTF-8}
\name{vegan-defunct}

\alias{commsimulator}
\alias{vegan-defunct}
%------ NOTE:  ../R/vegan-defunct.R   must be synchronized with this!
\title{Defunct Functions in Package \pkg{vegan}}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  These functions or variables listed here are no longer part of 
  \pkg{vegan}.
}
\usage{
commsimulator(x, method, thin=1)
}
\arguments{
  \item{x, method, thin}{Arguments to the function.}
}

\details{
  %-- Either:
  %There are currently no deprecated functions in this package.
  %-- Or:
  %% explain *why* it's deprecated, and \code{\link{..}} to new
  
  Function \code{commsimulator} is replaced with
  \code{\link{make.commsim}} which defines the Null models, and
  functions \code{\link{nullmodel}} and
  \code{\link{simulate.nullmodel}} that check the input data and
  generate the Null model communities.

}

\seealso{
  \code{\link{Defunct}}
}
\keyword{misc}
