\name{influence.cca}
\alias{hatvalues.cca}
\alias{hatvalues.rda}
\alias{sigma.cca}
\alias{rstandard.cca}
\alias{rstudent.cca}
\alias{cooks.distance.cca}
\alias{SSD.cca}
\alias{vcov.cca}
\alias{qr.cca}
\alias{df.residual.cca}

\title{Linear Model Diagnostics for Constrained Ordination}

\description{
  Influence statistics and other multiple linear model statistics for
  constrained ordination methods.
}

\usage{

\method{hatvalues}{cca}(model, ...)
\method{rstandard}{cca}(model, type = c("response", "canoco"), ...)
\method{rstudent}{cca}(model, type = c("response", "canoco"), ...)
\method{cooks.distance}{cca}(model, type = c("response", "canoco"), ...)

\method{sigma}{cca}(object, type = c("response", "canoco"), ...)
\method{vcov}{cca}(object, type = c("response", "canoco"), ...)
\method{SSD}{cca}(object, type = c("response", "canoco"), ...)

\method{qr}{cca}(x, ...)
\method{df.residual}{cca}(x, ...)

}

\arguments{

  \item{model, object, x}{A constrained ordination result object.}

  \item{type}{Type of statistics used for extracting raw residuals and
    residual standard deviation (\code{sigma}). Either
    \code{"response"} for species data or difference of WA and LC
    scores for \code{"canoco"}.}

}

\details{
  Linear model statistics (\code{\link{lm}}) methods for constrained
  ordination objects from \code{\link{cca}}, \code{\link{rda}},
  \code{\link{dbrda}}, and \code{\link{capscale}}.

}

\author{Jari Oksanen}

\seealso{Corresponding \code{\link{lm}} methods and
\code{\link{as.mlm.cca}}.}

\keyword{ models }
\keyword{ multivariate }