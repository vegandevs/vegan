\name{nobs.adonis}
\alias{nobs.adonis}
\alias{nobs.betadisper}
\alias{nobs.cca}
\alias{nobs.CCorA}
\alias{nobs.decorana}
\alias{nobs.isomap}
\alias{nobs.metaMDS}
\alias{nobs.pcnm}
\alias{nobs.procrustes}
\alias{nobs.rad}
\alias{nobs.varpart}
\alias{nobs.wcmdscale}

\title{
  Extract the Number of Observations from a vegan Fit.
}
\description{
  Extract the number of \sQuote{observations} from a \pkg{vegan} model fit.
}
\usage{
\method{nobs}{adonis}(object, ...)
}

\arguments{
  \item{object}{
    A fitted model object.
}
  \item{\dots}{
   Further arguments to be passed to methods.
}
}

\details{ Function \code{nobs} is generic in \R, and
  \pkg{vegan} provides methods for objects from \code{\link{adonis}},
  \code{\link{betadisper}}, \code{\link{cca}} and other related
  methods, \code{\link{CCorA}}, \code{\link{decorana}},
  \code{\link{isomap}}, \code{\link{metaMDS}}, \code{\link{pcnm}},
  \code{\link{procrustes}}, \code{\link{radfit}},
  \code{\link{varpart}} and \code{\link{wcmdscale}}.  }

\value{ A single number, normally an integer, giving the number of
  observations.  }

\author{
  Jari Oksanen
}

\keyword{models}

