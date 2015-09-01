\name{as.mlm.cca}
\alias{as.mlm}
\alias{as.mlm.cca}
\alias{as.mlm.rda}
\title{Refit Constrained Ordination as a Multiple Response Linear Model }
\description{
  Functions refit results of constrained ordination (\code{\link{cca}},
  \code{\link{rda}}, \code{\link{capscale}}) as a multiple response
  linear model (\code{\link{lm}}). This allows finding influence
  statistics (\code{\link{influence.measures}}).  This also allows
  deriving several other statistics, but most of these are biased and
  misleading, since refitting ignores a major component of variation in
  constrained ordination.
}
\usage{
as.mlm(x)
}
\arguments{
  \item{x}{Constrained ordination result. }
}
\details{
  Popular algorithm for constrained ordination is based on iteration
  with regression where weighted averages of sites are used as dependent
  variables and constraints as independent variables.
  Statistics of linear regression
  are a natural by-product in this algorithm. Constrained ordination in
  \pkg{vegan} uses different algorithm, but to obtain linear regression
  statistics you can refit an ordination result as a multiple response
  linear model (\code{\link{lm}}). This regression ignores residual
  unconstrained variation in the data, and therefore estimates of
  standard error are strongly biased and much too low. You can get
  statistics like \eqn{t}-values of coefficients, but you should not use
  these because of this bias. Some useful information you can get with
  refitted models are statistics for detecting influential observations
  (\code{\link{influence.measures}} including
  \code{\link{cooks.distance}}, \code{\link{hatvalues}}). 
}
\value{
  Function returns an object of multiple response linear model of class
  \code{"mlm"}  documented with \code{\link{lm}}.
}

\author{ Jari Oksanen }
\note{ 
  You can use these functions to find \eqn{t}-values of coefficients
  using \code{summary.mlm}, but you should not do this because the
  method ignores unconstrained residual variation. You also can find
  several other statistics for (multiple response) linear models with
  similar bias. This bias is not a unique feature in \pkg{vegan}
  implementation, but also applies to implementations in other
  software.

  Some statistics of linear models can be found without using
  these functions: \code{\link{coef.cca}} gives the regression
  coefficients, \code{\link{spenvcor}} the species-environment
  correlation, \code{\link{intersetcor}} the interset correlation,
  \code{\link{vif.cca}} the variance inflation factors.
}
\seealso{\code{\link{cca}}, \code{\link{rda}}, \code{\link{capscale}},
  \code{\link{cca.object}}, \code{\link{lm}}, \code{\link{summary.mlm}},
  \code{\link{influence.measures}}. }
\examples{
data(varespec)
data(varechem)
mod <- cca(varespec ~ Al + P + K, data=varechem)
lmod <- as.mlm(mod)
## Coefficients
lmod
coef(mod)
## Influential observations
influence.measures(lmod)
plot(mod, type = "n")
points(mod, cex = 10*hatvalues(lmod), pch=16, xpd = TRUE)
text(mod, display = "bp", col = "blue") 
}

\keyword{ models }
\keyword{ multivariate }
