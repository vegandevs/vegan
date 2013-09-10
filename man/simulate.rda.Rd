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
\method{simulate}{rda}(object, nsim = 1, seed = NULL, indx = NULL,
    rank = "full", correlated = FALSE, ...)
}
\arguments{

  \item{object}{an object representing a fitted \code{\link{rda}},
    \code{\link{cca}} or \code{\link{capscale}} model.}

  \item{nsim}{number of response matrices to be simulated. Only one
    dissimilarity matrix is returned for \code{\link{capscale}}, and
    larger \code{nsim} is an error.}

  \item{seed}{an object specifying if and how the random number
    generator should be initialized (\sQuote{seeded}). See 
    \code{\link{simulate}} for details. }

  \item{indx}{Index of residuals added to the fitted values, such as
    produced by \code{\link[permute]{shuffleSet}} or
    \code{\link{sample}}.  The index can have duplicate entries so
    that bootstrapping is allowed. If \code{nsim} \eqn{>1}, the output
    should be compliant with \code{\link[permute]{shuffleSet}} with
    one line for each simulation.  If \code{nsim} is missing, the
    number of rows of \code{indx} is used to define the number of
    simulations, but if \code{nsim} is given, it should match number
    of rows in \code{indx}. If null, parametric simulation is used and
    Gaussian error is added to the fitted values.}

  \item{rank}{The rank of the constrained component: passed to
    \code{\link{predict.rda}} or \code{\link{predict.cca}}. }

  \item{correlated}{Are species regarded as correlated in parametric
    simulation or when \code{indx} is not given? If
    \code{correlated = TRUE}, multivariate Gaussian random error is
    generated, and if \code{FALSE}, Gaussian random error is generated
    separately for each species. The argument has no effect in
    \code{\link{capscale}} which has no information on species.}

  \item{\dots}{additional optional arguments (ignored). }
}

\details{ The implementation follows \code{"lm"} method of
  \code{\link{simulate}}, and adds Gaussian (Normal) error to the fitted
  values (\code{\link{fitted.rda}}) using function \code{\link{rnorm}}
  if \code{correlated = FALSE} or \code{\link[MASS]{mvrnorm}} if
  \code{correlated = TRUE}. The standard deviations (\code{\link{rnorm}})
  or covariance matrices for species (\code{\link[MASS]{mvrnorm}}) are
  estimated from the residuals after fitting the constraints.
  Alternatively, the function can take a permutation index that is used
  to add permuted residuals (unconstrained component) to the fitted
  values. Raw data are used in \code{\link{rda}}. Internal Chi-square
  transformed data are used in \code{\link{cca}} within the function,
  but the returned matrix is similar to the original input data. The
  simulation is performed on internal metric scaling data in
  \code{\link{capscale}}, but the function returns the Euclidean
  distances calculated from the simulated data.  The simulation uses
  only the real components, and the imaginary dimensions are ignored.  }

\value{ If \code{nsim = 1}, returns a matrix or dissimilarities (in
  \code{\link{capscale}}) with similar additional arguments on random
  number seed as \code{\link{simulate}}. If \code{nsim > 1}, returns a
  similar array as returned by \code{\link{simulate.nullmodel}} with
  similar attributes.  }

\author{Jari Oksanen}


\seealso{ \code{\link{simulate}} for the generic case and for
  \code{\link{lm}} objects, and \code{\link{simulate.nullmodel}} for
  community null model simulation. Functions \code{\link{fitted.rda}}
  and \code{\link{fitted.cca}} return fitted values without the error
  component. See \code{\link{rnorm}} and \code{\link[MASS]{mvrnorm}}
  (\pkg{MASS} package) for simulating Gaussian random error. }

\examples{
data(dune)
data(dune.env)
mod <- rda(dune ~  Moisture + Management, dune.env)
## One simulation
update(mod, simulate(mod) ~  .)
## An impression of confidence regions of site scores
plot(mod, display="sites")
for (i in 1:5) lines(procrustes(mod, update(mod, simulate(mod) ~ .)), col="blue")
## Simulate a set of null communities with permutation of residuals
simulate(mod, indx = shuffleSet(nrow(dune), 99))
}
\keyword{ models }
\keyword{ datagen }
\keyword{ multivariate }
