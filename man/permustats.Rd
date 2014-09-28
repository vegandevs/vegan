\name{permustats}
\alias{permustats}
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
\alias{density.permustats}
\alias{qqnorm.permustats}
\alias{qqmath.permustats}

\title{
  Extract, Analyse and Display Permutation Results
}

\description{ 
  The \code{permustats} and \code{density} functions can directly
  access the permutation results of \pkg{vegan} functions.
}

\usage{
permustats(x, ...)
\method{summary}{permustats}(object, ...)
\method{densityplot}{permustats}(x, data, xlab = "Permutations", ...)
\method{density}{permustats}(x, observed = TRUE, ...)
\method{qqnorm}{permustats}(y, observed = TRUE, ...)
\method{qqmath}{permustats}(x, data, observed = TRUE, ylab = "Permutations", ...)
}

\arguments{
  \item{object, x, y}{The object to be handled.}
  
  \item{xlab, ylab}{Arguments of
    \code{\link[lattice]{densityplot}} and
    \code{\link[lattice]{qqmath}} functions.}
  \item{observed}{Add observed statistic among permutations.}
  \item{data}{Ignored.}
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

. However, you can use \code{densityplot} on
  \code{permustats} result for these methods.  NB, there are no
  \code{permustats} and \code{density} methods for
  \code{\link{anova.cca}}, but only for \code{\link{permutest.cca}}.

}

\value{
  The \code{permustats} function returns an object of class
  \code{"permustats"}. This is a list of items \code{"statistic"} for
  observed statistics, \code{permutations} which contains permuted
  values, and \code{alternative} which contains text defining the
  character of the test (\code{"two.sided"}, \code{"less"} or
  \code{"greater"}.).
}

\author{
  Jari Oksanen
}

\seealso{
  \code{\link{density}}, \code{\link[lattice]{densityplot}},
  \code{\link{qqnorm}}, \code{\link[lattice]{qqmath}}.
}

\examples{
data(dune)
data(dune.env)
mod <- adonis(dune ~ Management, data = dune.env)
## use permustats
perm <- permustats(mod)
summary(perm)
densityplot(perm)
qqmath(perm)
}

\keyword{ distribution }
\keyword{ smooth }
