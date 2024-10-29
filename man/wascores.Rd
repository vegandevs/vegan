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
\method{scores}{wascores}(x, display = c("wa", "stdev", "var", "se", "n2", "raw"),
             p = 0.95, ...)
}

\arguments{
  \item{x}{Environmental variables or ordination scores, or for
    \code{wascores} object with \code{stdev = TRUE}.}
  \item{w}{Weights: species abundances.}
  \item{expand}{Expand weighted averages so that they have the same
    weighted variance as the corresponding environmental variables.}
  \item{stdev}{Estimate weighted standard deviation of WA scores.}
  \item{display}{Type of scores returned.}
  \item{..., p}{Ignored.}
  
}
\details{

  Function \code{wascores} computes weighted averages. Weighted averages
  \dQuote{shrink}: they cannot be more extreme than values used for
  calculating the averages. With \code{expand = TRUE}, the function
  \dQuote{deshrinks} the weighted averages by making their biased
  weighted variance equal to the biased weighted variance of the
  corresponding environmental variable. There are many other methods of
  deshrinking, and the one used in this function is the simplest
  possible one and fixes only the shrinking of variance but not bias,
  like the popular regression based deshrinking.  Function
  \code{eigengrad} returns the inverses of squared expansion factors or
  the attribute \code{shrinkage} of the \code{wascores} result for each
  environmental gradient.  This is equal to the constrained eigenvalue
  of \code{\link{cca}} when only this one gradient was used as a
  constraint, and describes the strength of the gradient.

  With \code{stdev = TRUE} the function estimates the unbiased weighted
  standard deviation of the WA estimates using \code{\link{cov.wt}}. For
  unbiased standard deviation the virtual number of observations is
  equal to inverse Simpson index of diversity also known as Hill number
  N2. The numeric results can be accessed with \code{scores}
  function.

  The method is closely linked to correspondence analysis
  (\code{\link{ca}}, \code{\link{cca}}). Repeated use of \code{wascores}
  will converge to unconstrained correspondence analysis
  (\code{\link{ca}}) which is the reason for calling it Reciprocal
  Averaging. Function \code{\link{calibrate.cca}} for one-constraint
  constrained correspondence analysis (\code{\link{cca}}) will return
  weighted averages with different deshrinking.  Function
  \code{\link{tolerance}} uses the same algebra for weighted standard
  deviation, but bases the variance on linear combination
  scores which are scaled constraining variables instead of the weighted
  averages of the sites like \code{wascores}.

}

\value{
  If \code{stdev = TRUE}, function returns an object of class
  \code{"wascores"} with items
  
  \item{wa}{A matrix of weighted averages with.  If \code{expand=TRUE},
    attribute \code{shrinkage} has the inverses of squared expansion
    factors or \code{\link{cca}} eigenvalues for the variable.}
  \item{stdev}{a matrix of weighted standard deviations}
  \item{n2}{effective sample sizes which are equal to inverse Simpson
    diversity or Hill number N2}

  If \code{stdev = FALSE} (default), only the plain matrix \code{wa} is
  returned.  Function \code{eigengrad} returns only the \code{shrinkage}
  attribute.  }

\author{ Jari Oksanen }

\seealso{ \code{\link{cca}}, \code{\link{calibrate.cca}},
  \code{\link{tolerance.cca}}. WA scores are used for species in
  \code{\link{metaMDS}}. }

\examples{
data(varespec)
data(varechem)
vare.dist <- vegdist(wisconsin(varespec))
vare.mds <- monoMDS(vare.dist)
vare.points <- postMDS(vare.mds$points, vare.dist)
vare.wa <- wascores(vare.points, varespec)
plot(scores(vare.points), pch="+", asp=1)
text(vare.wa, rownames(vare.wa), cex=0.8, col="blue")
## Omit rare species (frequency <= 4)
freq <- apply(varespec>0, 2, sum)
plot(scores(vare.points), pch="+", asp=1)
text(vare.wa[freq > 4,], rownames(vare.wa)[freq > 4],cex=0.8,col="blue")
## Works for environmental variables, too.
head(spwa <- wascores(varechem, varespec, expand = TRUE))
## use species WA for Humdepth to estimate site WA
wascores(spwa[, "Humdepth"], t(varespec), expand = TRUE)
## use cca to get the weighted average of Humdepth (different deshrinking)
calibrate(cca(varespec ~ Humdepth, varechem))
## And the strengths of these variables are:
eigengrad(varechem, varespec)
}
\keyword{ multivariate }
\keyword{ univar }



