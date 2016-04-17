\name{RsquareAdj}
\alias{RsquareAdj}
\alias{RsquareAdj.default}
\alias{RsquareAdj.rda}
\alias{RsquareAdj.cca}
\alias{RsquareAdj.lm}
\alias{RsquareAdj.glm}

\Rdversion{1.1}

\title{
Adjusted R-square
}
\description{
  The functions finds the adjusted R-square.
}
\usage{
\method{RsquareAdj}{default}(x, n, m, ...)
\method{RsquareAdj}{rda}(x, ...)
\method{RsquareAdj}{cca}(x, permutations = 1000, ...)
}

\arguments{

  \item{x}{ Unadjusted R-squared or an object from which the terms for
  evaluation or adjusted R-squared can be found.}
  
  \item{n, m}{Number of observations and number of degrees of freedom
  in the fitted model.}
  
  \item{permutations}{Number of permutations to use when computing the adjusted 
  R-squared for a cca. The permutations can be calculated in parallel by
  specifying the number of cores which is passed to \code{\link{permutest}}}

  \item{\dots}{ Other arguments (ignored) except in the case of cca in 
  which these arguments are passed to \code{\link{permutest}}.} 
}

\details{ The default method finds the adjusted
  R-squared from the unadjusted R-squared, number of observations, and
  number of degrees of freedom in the fitted model. The specific
  methods find this information from the fitted result
  object. There are specific methods for \code{\link{rda}},
  \code{\link{cca}}, \code{\link{lm}} and \code{\link{glm}}. Adjusted,
  or even unadjusted, R-squared may not be available in some cases,
  and then the functions will return \code{NA}. There is no adjusted
  in partial \code{\link{rda}}, and
  R-squared values are available only for \code{\link{gaussian}}
  models in \code{\link{glm}}.

  The adjusted, \eqn{R^2}{R-squared} of \code{cca} is computed using a 
  permutation approach developed by Peres-Neto et al. (2006). By default 1000
  permutations are used.
  
  The raw \eqn{R^2}{R-squared} of partial \code{rda} gives the
  proportion explained after removing the variation due to conditioning
  (partial) terms; Legendre et al. (2011) call this semi-partial
  \eqn{R^2}{R-squared}. The adjusted \eqn{R^2}{R-squared} is found as
  the difference of adjusted \eqn{R^2}{R-squared} values of joint effect
  of partial and constraining terms and partial term alone, and it is
  the same as the adjusted \eqn{R^2}{R-squared} of component 
  \code{[a] = X1|X2} in two-component variation partition in 
  \code{\link{varpart}}.
  }

\value{ The functions return a list of items \code{r.squared} and
\code{adj.r.squared}.  
}

\references{
  Legendre, P., Oksanen, J. and ter Braak, C.J.F. (2011). Testing the
  significance of canonical axes in redundancy analysis. 
  \emph{Methods in Ecology and Evolution} 2, 269--277.
  
  Peres-Neto, P., P. Legendre, S. Dray and D. Borcard. 2006. Variation
  partitioning of species data matrices: estimation and comparison of
  fractions. \emph{Ecology} 87, 2614--2625.
}

\seealso{
  \code{\link{varpart}} uses \code{RsquareAdj}.
}
\examples{
data(mite)
data(mite.env)
## rda
m <- rda(decostand(mite, "hell") ~  ., mite.env)
RsquareAdj(m)
## cca
m <- cca(decostand(mite, "hell") ~  ., mite.env)
RsquareAdj(m)
## default method
RsquareAdj(0.8, 20, 5)
}
\keyword{ univar }
\keyword{ multivariate }
