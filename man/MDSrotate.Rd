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

  \item{vec}{ A continuous environmental variable or a matrix of such
    variables. The number of variables must be lower than the number of
    dimensions, and the solution is rotated to these variables in the
    order they appear in the matrix.}

  \item{na.rm}{ Remove missing values from the continuous variable
    \code{vec}.}

  \item{\dots}{ Other arguments (ignored). }

}

\details{ The orientation and rotation are undefined in multidimensional
  scaling.  Functions \code{\link{metaMDS}} and \code{\link{metaMDS}}
  can rotate their solutions to principal components so that the
  dispersion of the points is highest on the first dimension. Sometimes
  a different rotation is more intuitive, and \code{MDSrotate} allows
  rotation of the result so that the first axis is parallel to a given
  external variable or two first variables are completely in a
  two-dimensional plane etc. If several external variables are supplied,
  they are applied in the order they are in the matrix. First axis is
  rotated to the first supplied variable, and the second axis to the
  second variable. Because variables are usually correlated, the second
  variable is not usually aligned with the second axis, but it is
  uncorrelated to later dimensions. There must be at least one free
  dimension: the number of external variables must be lower than the
  number of dimensions, and all used environmental variables are
  uncorrelated with that free dimension.}

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
ef <- envfit(mod ~ pH, varechem, permutations = 0)
plot(ef)
ordisurf(mod ~ pH, varechem, knots = 1, add = TRUE)
}

\keyword{multivariate}

