\name{permustats}
\alias{permustats.adonis}
\alias{permustats.anosim}
\alias{permustats.CCorA}
\alias{permustats.envfit}
\alias{permustats.factorfit}
\alias{permustats.mantel}
\alias{permustats.mrpp}
\alias{permustats.mso}
\alias{permustats.oecosimu}
\alias{permustats.permutest.betadisper}
\alias{permustats.permutest.cca}
\alias{permustats.protest}
\alias{permustats.vectorfit}
\alias{summary.permustats}
\alias{densityplot.permustats}
\alias{density.adonis}
\alias{density.anosim}
\alias{density.mantel}
\alias{density.mrpp}
\alias{density.permutest.cca}
\alias{density.protest}
\alias{plot.vegandensity}
\alias{densityplot.adonis}

\title{
  Extract, Analyse and Display Densities of Permutation Results
}

\description{ 
  The \code{permustats} and \code{density} functions can directly
  access the permutation results of \pkg{vegan} functions.
}

\usage{
permustats(x, ...)
\method{summary}{permustats}(object, ...)
\method{densityplot}{permustats}(x, xlab = "Permutations", ...)
\method{density}{adonis}(x, ...)
\method{plot}{vegandensity}(x, main = NULL, xlab = NULL, ylab = "Density", 
   type = "l", zero.line = TRUE, obs.line = TRUE, ...)
\method{densityplot}{adonis}(x, data, xlab = "Null", ...)
}

\arguments{
  \item{object, x}{The object to be handled.}
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

  The \code{permustats} function extracts permutation results and
  observed statistics.  The \code{summary} estimate the \eqn{z}
  values, also known as standardized effect size (SES) as the
  difference of obseved statistic and mean of permutations divided by
  the standard deviation of the permutations. It also prints the the
  mean, median, and limits within which are 95\% of permutations. For
  two-sided test these are (2.5\%, 97.5\%) and for one-sided tests
  either 5\% or 95\% quantile depending on the text direction. The
  \code{densityplot} method displays the density of permutations and
  marks the observed statistic as a vertical line (see
  \code{\link[lattice]{densityplot}}). The observed statistic is
  always put among the permuted values so that the plot is consistnet
  with significance tests. The \code{permustats} can extract
  permutation statistics from the results of \code{\link{adonis}},
  \code{\link{anosim}}, \code{\link{mantel}},
  \code{\link{mantel.partial}}, \code{\link{mrpp}},
  \code{\link{oecosimu}}, \code{\link{permutest.cca}}, and
  \code{\link{protest}}.

  The \code{density} function can directly access permutation results
  of the same function as \code{permustats}.  The \code{density}
  function is identical to \code{\link{density.default}} and takes all
  its arguments, but adds the observed statistic to the result as item
  \code{"observed"}. The observed statistic is also put among the
  permuted values so that the results are consistent with significance
  tests. The \code{plot} method is similar to the default
  \code{\link{plot.density}}, but can also add the observed statistic
  to the graph as a vertical line.  In \code{\link{adonis}} it is also
  possible to use direclty \code{densityplot} function.

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

  Functions that can return several permuted statistics simultaneously
  also have \code{\link[lattice]{densityplot}} method
  (\code{\link{adonis}}, \code{\link{oecosimu}} and diversity
  partitioning functions based on \code{oecosimu}).  The standard
  \code{\link{density}} can only handle univariate data, and a warning
  is issued if the function is used for a model with several observed
  statistics. However, you can use \code{densityplot} on
  \code{permustats} result for these methods.  NB, there are no
  \code{permustats} and \code{density} methods for
  \code{\link{anova.cca}}, but only for \code{\link{permutest.cca}}.

}

\value{
  The \code{permustats} function return an object of class
  \code{"permustats"}. This is a list of items \code{"statistic"} for
  observed statistics, \code{permutations} which contains permuted
  values, and \code{alternative} which contains text defining the
  character of the test (\code{"two.sided"}, \code{"less"} or
  \code{"greater"}.).

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
## use permustats
perm <- permustats(mod)
summary(perm)
densityplot(perm)
## direct use of density
plot(density(mod))
mod <- adonis(dune ~ Management * Moisture, dune.env)
densityplot(mod)
}

\keyword{ distribution }
\keyword{ smooth }
