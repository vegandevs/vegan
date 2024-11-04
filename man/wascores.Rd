\name{wascores}
\alias{wascores}
\alias{scores.wascores}
\alias{eigengrad}

\title{ Weighted Averages Scores for Species }
\description{
  Computes Weighted Averages scores of species for ordination
  configuration or for environmental variables.
}
\usage{
wascores(x, w, expand = FALSE, stdev = FALSE)
eigengrad(x, w)
\method{scores}{wascores}(x, display = c("wa", "stdev", "var", "se", "n2", "raw"), ...)
}

\arguments{
  \item{x}{Environmental variables or ordination scores, or for
    \code{wascores} object with \code{stdev = TRUE}.}
  \item{w}{Weights: species abundances.}
  \item{expand}{Expand weighted averages so that they have the same
    weighted variance as the corresponding environmental variables.}
  \item{stdev}{Estimate weighted standard deviation of WA scores.}
  \item{display}{Type of scores returned.}
  \item{\dots}{Other arguments passed to functions (currently ignored).}
  
}
\details{

  Weighted Averages are a classical way of estimating the species
  optima along continuous environmental variables
  (a.k.a. gradients). Function \code{wascores} is a simple function
  that is mainly designed to add species scores to unimodal
  ordinations (\code{\link{metaMDS}}, \code{\link{sppscores}}) or
  ordering rows or columns to give diagonal pattern of tabulation
  (\code{\link{vegemite}}, \code{\link{tabasco}}). It can also be used
  to find species \dQuote{optima} or sampling unit calibrations for
  community data. For this purpose, specialized packages such
  \pkg{analogue} are recommended (but see
  \code{\link{calibrate.cca}}).

  First argument of \code{wascores} is the variable or a matrix of
  variables for which weighted averages are needed, and the second
  argument is the matrix of weights. In classical approaches weights are
  a community matrix, where taxon abundances define the weights. The
  number of rows must match. If the first argument is for taxa
  (columns), community weight matrix must be transposed. 

  Weighted averages \dQuote{shrink}: they cannot be more extreme than
  values used for calculating the averages. With \code{expand = TRUE},
  the function \dQuote{deshrinks} the weighted averages making their
  weighted variance equal to the weighted variance of the
  corresponding input variable. Specialized packages (such as
  \pkg{analogue}) offer a wider range of deshrinking alternatives, but
  deshrinking can also made after the analysis (see
  Examples). Function \code{eigengrad} returns the strength of
  expansion as attribute \code{shrinkage} of the \code{wascores}
  result for each environmental gradient.  The shrinkage equal to the
  constrained eigenvalue of \code{\link{cca}} when only this one
  gradient was used as a constraint, and describes the strength of the
  gradient.

  With \code{stdev = TRUE} the function estimates the unbiased
  weighted standard deviation of the WA estimates using
  \code{\link{cov.wt}}. For unbiased standard deviation the virtual
  number of observations is equal to inverse Simpson index of
  diversity also known as Hill number N2 (see
  \code{\link{diversity}}). The numeric results can be accessed with
  \code{scores} function.  Function \code{\link{tolerance}} uses the
  same algebra for weighted standard deviation, but bases the variance
  on linear combination scores (constaints) variables instead of the
  weighted averages of the sites like \code{wascores}.

  Weighted averages are closely linked to correspondence analysis
  (\code{\link{ca}}, \code{\link{cca}}). Repeated use of
  \code{wascores} will converge to the first axis of unconstrained
  correspondence analysis (\code{\link{ca}}) which therefore is also
  known as Reciprocal Averaging (Hill 1973). Constrained
  correspondence analysis (\code{\link{cca}}) is equivalent to
  weighted averages and \code{\link{calibrate.cca}} will return
  weighted averages of the constraint with different deshrinking.

}

\value{
  If \code{stdev = TRUE}, function returns an object of class
  \code{"wascores"} with items
  
  \item{wa}{A matrix of weighted averages with.  If
    \code{expand=TRUE}, attribute \code{shrinkage} has the inverses of
    squared expansion factors or \code{\link{cca}} eigenvalues for the
    variable and attribute \code{centre} for the weighted means of the
    variables.}

  \item{stdev}{a matrix of weighted standard deviations}
  \item{n2}{effective sample sizes which are equal to inverse Simpson
    diversity or Hill number N2}

  If \code{stdev = FALSE} (default), only the plain matrix \code{wa}
  is returned.  Function \code{eigengrad} returns only the
  \code{shrinkage} attribute. With \code{stdev = TRUE} only a brief
  summary of the result is printed, and the individvual scores can be
  accessed with \code{scores} function. }

\author{ Jari Oksanen }

\seealso{\code{\link{calibrate.cca}}, \code{\link{tolerance.cca}},
  \code{\link{sppscores}}.}

\references{

  Hill, M.O. (1973) Reciprocal averaging: An eigenvector method of
  ordination. \emph{Journal of Ecology} 61, 237--249.
  
}

\examples{
data(mite, mite.env)
## add species points to ordination
mod <- monoMDS(vegdist(mite))
plot(mod)
## add species points; sppscores does the same and can also add the
## species scores to mod
points(wascores(scores(mod), mite, expand = TRUE), pch="+", col=2)
## Get taxon optima for WatrCont
head(wascores(mite.env$WatrCont, mite))
## WA calibration: site WA from species WA; NB using transpose for site WA
spwa <- wascores(mite.env$WatrCont, mite, expand = TRUE)
wacalib <- wascores(spwa, t(mite), expand = TRUE)
plot(wacalib ~ WatrCont, data=mite.env)
abline(0, 1)
## use traditional 'inverse' regression deshrinking instead of wascores
## 'expand'
wareg <- fitted(lm(WatrCont ~ wacalib, data=mite.env))
head(cbind("WatrCont" = mite.env$WatrCont, "expand" = drop(wacalib),
    "regression" = wareg))
## Reciprocal Averaging algorithm for Correspondence Analysis
## start with random values
u <- runif(nrow(mite))
## repeat the following steps so long that the shrinkage converges
v <- wascores(u, mite, expand = TRUE)
u <- wascores(v, t(mite), expand = TRUE)
attr(u, "shrinkage") # current estimate of eigenvalue
## The strengths of two continuous variables in the data set
eigengrad(mite.env[, 1:2], mite)
}
\keyword{ multivariate }
\keyword{ univar }



