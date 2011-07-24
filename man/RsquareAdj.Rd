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
}

\arguments{

  \item{x}{ Unadjusted R-squared or an object from which the terms for
  evaluation or adjusted R-squared can be found.}
  
  \item{n, m}{Number of observations and number of degrees of freedom
  in the fitted model.}

  \item{\dots}{ Other arguments (ignored).} 
}

\details{ The default method finds the adjusted
  R-squared from the unadjusted R-squared, number of observations, and
  number of degrees of freedom in the fitted model. The specific
  methods find this information from the fitted result
  object. There are specific methods for \code{\link{rda}},
  \code{\link{cca}}, \code{\link{lm}} and \code{\link{glm}}. Adjusted,
  or even unadjusted, R-squared may not be available in some cases,
  and then the functions will return \code{NA}. There is no adjusted
  R-squared in \code{\link{cca}}, in partial \code{\link{rda}}, and
  R-squared values are available only for \code{\link{gaussian}}
  models in \code{\link{glm}}. 
}

\value{ The functions return a list of items \code{r.squared} and
\code{adj.r.squared}.  
}

\references{ 
  Peres-Neto, P., P. Legendre, S. Dray and D. Borcard. 2006. Variation
  partitioning of species data matrices: estimation and comparison of
  fractions. \emph{Ecology} 87: 2614-2625.  }

\seealso{
  \code{\link{varpart}} uses \code{RsquareAdj}.
}
\examples{
data(mite)
data(mite.env)
## rda
m <- rda(decostand(mite, "hell") ~  ., mite.env)
RsquareAdj(m)
## default method
RsquareAdj(0.8, 20, 5)
}
\keyword{ univar }
\keyword{ multivariate }
