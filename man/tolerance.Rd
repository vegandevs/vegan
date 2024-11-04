\name{tolerance}
\alias{tolerance}
\alias{tolerance.cca}
\alias{tolerance.decorana}

\title{Species tolerances and sample heterogeneities}
\usage{
tolerance(x, \dots)

\method{tolerance}{cca}(x, choices = 1:2, which = c("species","sites"),
          scaling = "species", useN2 = TRUE, hill = FALSE, \dots)

\method{tolerance}{decorana}(x, data, choices = 1:4,
          which = c("sites", "species"), useN2 = TRUE, ...)

}
\description{
  Species tolerances and sample heterogeneities.
}
\details{
  Function to compute species tolerances and site heterogeneity measures
  from unimodal ordinations (CCA & CA). Implements Eq 6.47 and 6.48 from
  the Canoco 4.5 Reference Manual (pages 178--179).

  Function \code{wascores} with \code{stdev = TRUE} uses the same
  algebra, but bases the standard deviations on weighted averages
  scores instead of linear combinations scores of \code{tolerance}.
}
\value{
  Matrix of tolerances/heterogeneities with some additional
  attributes: \code{which}, \code{scaling}, and \code{N2}, the latter of
  which will be \code{NA} if \code{useN2 = FALSE} or \code{N2} could not
  be estimated.
}
\author{Gavin L. Simpson and Jari Oksanen (\code{decorana} method).}
\arguments{
  \item{x}{object of class \code{"cca"}.}
  \item{choices}{numeric; which ordination axes to compute
    tolerances and heterogeneities for. Defaults to axes 1 and 2.}
  \item{which}{character; one of \code{"species"} or \code{"sites"},
    indicating whether species tolerances or sample heterogeneities
    respectively are computed.}
  \item{scaling}{character or numeric; the ordination scaling to
    use. See \code{\link{scores.cca}} for details.}
  \item{hill}{logical; if \code{scaling} is a character,
    these control whether Hill's scaling is used for (C)CA
    respectively. See \code{\link{scores.cca}} for details.}
  \item{useN2}{logical; should the bias in the tolerances /
    heterogeneities be reduced via scaling by Hill's N2?}
  \item{data}{Original input data used in \code{\link{decorana}}.  If
    missing, the function tries to get the same data as used in
    \code{decorana} call.}
  \item{\dots}{arguments passed to other methods.}
}
\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ ., data = dune.env)

## defaults to species tolerances
tolerance(mod)

## sample heterogeneities for CCA axes 1:6
tolerance(mod, which = "sites", choices = 1:6)
## average should be 1 with scaling = "sites", hill = TRUE
tol <- tolerance(mod, which = "sites", scaling = "sites", hill = TRUE,
   choices = 1:4)
colMeans(tol)
apply(tol, 2, sd)
## Rescaling tries to set all tolerances to 1
tol <- tolerance(decorana(dune))
colMeans(tol)
apply(tol, 2, sd)
}
