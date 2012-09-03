\name{density.adonis}
\alias{density.adonis}
\alias{density.anosim}
\alias{density.mantel}
\alias{density.mrpp}
\alias{density.permutest.cca}
\alias{density.protest}
\alias{plot.vegandensity}
\alias{densityplot.adonis}

\title{
  Kernel Density Estimation for Permutation Results in Vegan
}

\description{ 
  The \code{density} functions can directly access the permutation
  results of \pkg{vegan} functions, and \code{plot} can display the
  densities. The \code{densityplot} method can access and display the
  permutation results of functions that return permutations of several
  statistics simultaneously.  
}

\usage{
\method{density}{adonis}(x, ...)
\method{plot}{vegandensity}(x, main = NULL, xlab = NULL, ylab = "Density", 
   type = "l", zero.line = TRUE, obs.line = TRUE, ...)
}

\arguments{
  \item{x}{The object to be handled. For \code{density} and
     \code{densityplot} this is an object containing permutations. For
     \code{plot} this is a result of \pkg{vegan} \code{density}
     function.}
  \item{main, xlab, ylab, type, zero.line}{Arguments of
    \code{\link{plot.density}} and \code{\link[lattice]{densityplot}}
    functions.}
  \item{obs.line}{Draw vertical line for the observed
    statistic. Logical value \code{TRUE} draws a red line, and
    \code{FALSE} draws nothing. Alternatively, \code{obs.line} can be a
    definition of the colour used for the line, either as a numerical
    value from the \code{\link[grDevices]{palette}} or as the name of
    the colour, or other normal definition of the colour.}
  \item{\dots}{ Other arguments passed to the function. In
    \code{density} these are passed to \code{\link{density.default}}.}
}

\details{ 

  The \code{density} and \code{densityplot} function can directly access
  permutation results of most \pkg{vegan} functions.  The \code{density}
  function is identical to \code{\link{density.default}} and takes all
  its arguments, but adds the observed statistic to the result as item
  \code{"observed"}. The observed statistic is also put among the
  permuted values so that the results are consistent with significance
  tests. The \code{plot} method is similar to the default
  \code{\link{plot.density}}, but can also add the observed statistic to
  the graph as a vertical line.  The \code{densityplot} function is
  based on the same function in the \pkg{lattice} package (see
  \code{\link[lattice]{densityplot}}).

  The density methods are available for \pkg{vegan} functions
  \code{\link{adonis}}, \code{\link{anosim}}, \code{\link{mantel}},
  \code{\link{mantel.partial}}, \code{\link{mrpp}},
  \code{\link{permutest.cca}}, and \code{\link{protest}}.  The
  \code{density} function for \code{\link{oecosimu}} is documented
  separately, and it is also used for \code{\link{adipart}},
  \code{\link{hiersimu}} and \code{\link{multipart}}.

  All \pkg{vegan} \code{density} functions return an object of class
  \code{"vegandensity"} inheriting from \code{\link{density}}, and can
  be plotted with its \code{plot} method.  This is identical to the
  standard \code{plot} of \code{densiy} objects, but can also add a
  vertical line for the observed statistic.

  Functions that can return several permuted statistics simultaenously
  also have \code{\link[lattice]{densityplot}} method
  (\code{\link{adonis}}, \code{\link{oecosimu}} and diversity partioning
  functions based on \code{oecosimu}).  The standard
  \code{\link{density}} can only handle univariate data, and a warning
  is issued if the function is used for a model with several observed
  statistics.  The \code{\link[lattice]{densityplot}} method is available
  for \code{\link{adonis}} and \code{\link{oecosimu}} (documented
  separately). NB, there is no \code{density} method for
  \code{\link{anova.cca}}, but only for \code{\link{permutest.cca}}.

}

\value{
  The \code{density} function returns the standard \code{\link{density}}
  result object with one new item: \code{"observed"} for the observed
  value of the statistic. The functions have a specific \code{plot}
  method, but otherwise they use methods for
  \code{\link{density.default}}, such as \code{print} and \code{lines}.
}

\author{
  Jari Oksanen
}

\seealso{
  \code{\link{density.default}}.
}

\examples{
data(dune)
data(dune.env)
mod <- adonis(dune ~ Management, data = dune.env)
plot(density(mod))
library(lattice)
mod <- adonis(dune ~ Management * Moisture, dune.env)
densityplot(mod)
}

\keyword{ distribution }
\keyword{ smooth }
