\name{simulate.rda}
\alias{simulate.rda}
\alias{simulate.cca}
\alias{simulate.capscale}
\title{ Simulate Responses with Gaussian Error or Permuted Residuals for Constrained Ordination }

\description{ Function simulates a response data frame so that it adds
 Gaussian error to the fitted responses of Redundancy Analysis
 (\code{\link{rda}}), Constrained Correspondence Analysis
 (\code{\link{cca}}) or distance-based RDA (\code{\link{capscale}}).
 The function is a special case of generic \code{\link{simulate}}, and
 works similarly as \code{simulate.lm}.  }

\usage{
\method{simulate}{rda}(object, nsim = 1, seed = NULL, indx = NULL, rank = "full", ...)
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
    \code{\link{shuffle}} or \code{\link{sample}}. The index can
    have duplicate entries so that bootstrapping is allowed. If null,
    parametric simulation is used and Gaussian error is added to the
    fitted values.}
  \item{rank}{The rank of the constrained component: passed to
    \code{\link{predict.rda}} or \code{\link{predict.cca}}. }
  \item{\dots}{additional optional arguments (ignored). }
}

\details{ The implementation follows \code{"lm"} method of
  \code{\link{simulate}}, and adds Gaussian (Normal) error to the
  fitted values (\code{\link{fitted.rda}}) using function
  \code{\link{rnorm}}. The standard deviations are estimated
  independently for each species (column) from the residuals after
  fitting the constraints. Alternatively, the function can take a
  permutation index that is used to add permuted residuals
  (unconstrained component) to the fitted values. Raw data are used in
  \code{\link{rda}}. Internal Chi-square transformed data in
  \code{\link{cca}} within the function, but the returned data frame is 
  similar to the original input data. The simulation is performed on
  internal metric scaling data in \code{\link{capscale}}, but the
  function returns the Euclidean distances calculated from the simulated
  data.  The simulation uses only the real components, and the imaginary
  dimensions are ignored.
}

\value{ Returns a data frame with similar additional arguments on
  random number seed as \code{\link{simulate}}.  }

\author{Jari Oksanen}


\seealso{ \code{\link{simulate}} for the generic case and for
  \code{\link{lm}} objects. Functions \code{\link{fitted.rda}} and
  \code{\link{fitted.cca}} return fitted values without the error
  component.  
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
