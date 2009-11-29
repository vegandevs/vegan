\name{simulate.rda}
\alias{simulate.rda}
\alias{simulate.cca}
\alias{simulate.capscale}
\title{ Simulate Responses with Gaussian Error for Redundancy Analysis }

\description{ Function simulates a response data frame so that it adds
 Gaussian error to the fitted responses of Redundancy Analysis
 (\code{\link{rda}}). The function is a special case of generic
 \code{\link{simulate}}, and works similarly as \code{simulate.lm}. 
}

\usage{
\method{simulate}{rda}(object, nsim = 1, seed = NULL, indx = NULL, ...)
}
\arguments{
  \item{object}{an object representing a fitted \code{\link{rda}} model.}
  \item{nsim}{number of response vectors to simulate. (Not yet used, and 
    values above 1 will give an error). }
  \item{seed}{an object specifying if and how the random number
    generator should be initialized (\sQuote{seeded}). See 
    \code{\link{simulate}} for details. }
  \item{indx}{Index of residuals added to the fitted values, such as
    produced by  \code{\link{permuted.index}},
    \code{\link{permuted.index2}} or \code{\link{sample}}. The index can
    have duplicate entries so that bootstrapping is allowed. If null,
    parametric simulation is used and Gaussian error is added to the
    fitted values.}
  \item{\dots}{additional optional arguments (ignored). }
}

\details{ The implementation follows \code{"lm"} method of
  \code{\link{simulate}}, and adds Gaussian (Normal) error to the
  fitted values (\code{\link{fitted.rda}} using function
  \code{\link{rnorm}}. The standard deviations are estimated
  independently for each species (column) from the residuals after
  fitting the constraints.
}

\value{ Returns a data frame with similar additional arguments on
  random number seed as \code{\link{simulate}}.  }

\author{Jari Oksanen}

\note{ The function is not implemented for \code{\link{cca}} or
  \code{\link{capscale}} objects, but only for \code{\link{rda}}.  
}

\seealso{ \code{\link{simulate}} for the generic case and for
  \code{\link{lm}} objects. Function \code{\link{fitted.rda}} returns
  fitted values without the error component. 
}

\examples{
data(dune)
data(dune.env)
mod <- rda(dune ~  Moisture + Management, dune.env)
## One simulation
update(mod, simulate(mod) ~  .)
## An impression of confidence regions of site scores
plot(mod, display="sites")
for (i in 1:5) lines(procrustes(mod, update(mod, simulate(mod) ~ .)), col="blue")
}
\keyword{ models }
\keyword{ datagen }
\keyword{ multivariate }
