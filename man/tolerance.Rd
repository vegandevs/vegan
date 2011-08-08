\name{tolerance}
\alias{tolerance}
\alias{tolerance.cca}
\alias{print.tolerance.cca}
\title{Species tolerances and sample heterogeneities}
\usage{
tolerance(x, \dots)

\method{tolerance}{cca}(x, choices = 1:2, which = c("species","sites"),
    scaling = 2, useN2 = FALSE, \dots)
}
\description{
  Species tolerances and sample heterogeneities.
}
\details{
  Function to compute species tolerances and site heterogeneity measures
  from unimodal ordinations (CCA & CA). Implements Eq 6.47 and 6.48 from
  the Canoco 4.5 Reference Manual (pages 178-179).
}
\value{
  Matrix of tolerances/heterogeneities with some additional
  attributes.
}
\author{Gavin L. Simpson}
\arguments{
  \item{x}{object of class \code{"cca"}.}
  \item{choices}{numeric; which ordination axes to compute
    tolerances and heterogeneities for. Defaults to axes 1 and 2.}
  \item{which}{character; one of \code{"species"} or \code{"sites"},
    indicating whether species tolerances or sample heterogeneities
    respectively are computed.}
  \item{scaling}{numeric; the ordination scaling to use.}
  \item{useN2}{logical; should the bias in the tolerances /
    heterogeneities be reduced via scaling by Hill's N2?}
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
}
