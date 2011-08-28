\name{MDSrotate}
\alias{MDSrotate}

\title{
  Rotate First MDS Dimension Parallel to an External Variable
}

\description{ Function rotates a multidimensional scaling result so
  that its first dimension is parallel to an external (environmental
  variable). The function can handle the results from
  \code{\link{metaMDS}} or \code{\link{monoMDS}} functions.  }

\usage{
MDSrotate(object, vec, na.rm = FALSE, ...)
}

\arguments{
 
 \item{object}{ A result object from \code{\link{metaMDS}} or
    \code{\link{monoMDS}}.}

  \item{vec}{ A continuous environmental variable (vector of the same
    length as the number of points).}

  \item{na.rm}{ Remove missing values from the continuous variable
    \code{vec}.}

  \item{\dots}{ Other arguments (ignored). }

}

\details{ The orientation and rotation are undefined in
  multidimensional scaling.  Functions \code{\link{metaMDS}} and
  \code{\link{metaMDS}} can rotate their solutions to principal
  components so that the dispersion of the points is highest on the
  first dimension. Sometimes a different rotation is more intuitive,
  and \code{MDSrotate} allows rotation of the result so that the first
  axis is parallel to a given external variable.  
}

\value{ Function returns the original ordination result, but with
  rotated scores (both site and species if available), and the
  \code{pc} attribute of scores set to \code{FALSE}.  
}


\author{
  Jari Oksanen
}

\seealso{
  \code{\link{metaMDS}}, \code{\link{monoMDS}}.
}
\examples{
data(varespec)
data(varechem)
mod <- monoMDS(vegdist(varespec))
mod <- with(varechem, MDSrotate(mod, pH))
plot(mod)
ef <- envfit(mod ~ pH, varechem)
plot(ef)
ordisurf(mod ~ pH, varechem, knots = 1, add = TRUE)
}

\keyword{multivariate}

