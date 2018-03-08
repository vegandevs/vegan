\name{wascores}
\alias{wascores}
\alias{eigengrad}

\title{ Weighted Averages Scores for Species }
\description{
  Computes Weighted Averages scores of species for ordination
  configuration or for environmental variables.
}
\usage{
wascores(x, w, expand=FALSE)
eigengrad(x, w)
}

\arguments{
  \item{x}{Environmental variables or ordination scores.}
  \item{w}{Weights: species abundances.}
  \item{expand}{Expand weighted averages so that they have the same
    weighted variance as the corresponding environmental variables.  }
}
\details{
  Function \code{wascores} computes weighted averages. Weighted averages
  \dQuote{shrink}: they cannot be more extreme than values used for
  calculating the averages. With \code{expand = TRUE}, the function
  \dQuote{deshrinks} the weighted averages by making their biased
  weighted variance equal to the biased weighted variance of the
  corresponding environmental variable.  Function \code{eigengrad}
  returns the inverses of squared expansion factors or the attribute
  \code{shrinkage} of the \code{wascores} result for each environmental
  gradient.  This is equal to the constrained eigenvalue of
  \code{\link{cca}} when only this one gradient was used as a
  constraint, and describes the strength of the gradient.
}
\value{
  Function \code{wascores} returns a matrix where species define rows
  and ordination axes or environmental variables define columns. If
  \code{expand = TRUE}, attribute \code{shrinkage} has the inverses of
  squared expansion factors or \code{\link{cca}} eigenvalues for the
  variable.  Function \code{eigengrad} returns only the \code{shrinkage}
  attribute.
}

\author{ Jari Oksanen }

\seealso{ \code{\link{monoMDS}}, \code{\link{cca}}. }

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
wascores(varechem, varespec)
## And the strengths of these variables are:
eigengrad(varechem, varespec)
}
\keyword{ multivariate }
\keyword{ univar }



