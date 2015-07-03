\encoding{UTF-8}
\name{permustats}
\alias{permustats}
\alias{permustats.adonis}
\alias{permustats.anosim}
\alias{permustats.anova.cca}
\alias{permustats.CCorA}
\alias{permustats.envfit}
\alias{permustats.factorfit}
\alias{permustats.mantel}
\alias{permustats.mrpp}
\alias{permustats.mso}
\alias{permustats.oecosimu}
\alias{permustats.ordiareatest}
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
  The \code{permustats} function extracts permutation results of
  \pkg{vegan} functions. Its support functions can find quantiles and
  standardized effect sizes, plot densities and Q-Q plots.
}

\usage{
permustats(x, ...)
\method{summary}{permustats}(object, interval = 0.95, ...)
\method{densityplot}{permustats}(x, data, xlab = "Permutations", ...)
\method{density}{permustats}(x, observed = TRUE, ...)
\method{qqnorm}{permustats}(y, observed = TRUE, ...)
\method{qqmath}{permustats}(x, data, observed = TRUE, ylab = "Permutations", ...)
}

\arguments{
  \item{object, x, y}{The object to be handled.}
  \item{interval}{numeric; the coverage interval reported.}
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
  observed statistics from several \pkg{vegan} functions that perform
  permutations or simulations.  

  The \code{summary} method of \code{permustats} estimates the
  standardized effect sizes (SES) as the difference of observed
  statistic and mean of permutations divided by the standard deviation
  of permutations (also known as \eqn{z}-values). It also prints the
  the mean, median, and limits which contain \code{interval} percent
  of permuted values. With the default (\code{interval = 0.95}), for
  two-sided test these are (2.5\%, 97.5\%) and for one-sided tests
  either 5\% or 95\% quantile depending on the test direction. The
  mean, quantiles and \eqn{z} values are evaluated from permuted
  values without observed statistic.

  The \code{density} and \code{densityplot} methods display the
  kernel density estimates of permuted values. When observed value of
  the statistic is included in the permuted values, the
  \code{densityplot} method marks the observed statistic as a vertical
  line. However the \code{density} method uses its standard \code{plot}
  method and cannot mark the obseved value.

  The \code{qqnorm} and \code{qqmath} display Q-Q plots of
  permutations, optionally together with the observed value (default)
  which is shown as horizontal line in plots. \code{qqnorm} plots
  permutation values against standard Normal variate. \code{qqmath}
  defaults to the standard Normal as well, but can accept other
  alternatives (see standard \code{\link[lattice]{qqmath}}).

  Functions \code{\link{density}} and \code{\link{qqnorm}} are based on
  standard \R methods and accept their arguments. They only handle one
  statistic, and cannot be used when several test statistic were
  evaluated. The \code{\link[lattice]{densityplot}} and
  \code{\link[lattice]{qqmath}} are \pkg{lattice} graphics, and can be
  used both for one and several statistics.  All these functions pass
  arguments to their underlying functions; see their documentation.

  The \code{permustats} can extract permutation statistics from the
  results of \code{\link{adonis}}, \code{\link{anosim}}, 
  \code{\link{anova.cca}}, \code{\link{mantel}}, \code{\link{mantel.partial}},
  \code{\link{mrpp}}, \code{\link{oecosimu}}, \code{\link{ordiareatest}},
  \code{\link{permutest.cca}}, \code{\link{protest}}, and
  \code{\link{permutest.betadisper}}. 

}

\value{
  The \code{permustats} function returns an object of class
  \code{"permustats"}. This is a list of items \code{"statistic"} for
  observed statistics, \code{permutations} which contains permuted
  values, and \code{alternative} which contains text defining the
  character of the test (\code{"two.sided"}, \code{"less"} or
  \code{"greater"}). The \code{\link{qqnorm}} and
  \code{\link{density}} methods return their standard result objects.
}

\author{
  Jari Oksanen with contributions from Gavin L. Simpson
  (\code{permustats.permutest.betadisper} method and related
  modifications to \code{summary.permustats} and the \code{print}
  method) and Eduard \enc{Szöcs}{Szoecs} (\code{permustats.anova.cca).}
}

\seealso{
  \code{\link{density}}, \code{\link[lattice]{densityplot}},
  \code{\link{qqnorm}}, \code{\link[lattice]{qqmath}}.
}

\examples{
data(dune)
data(dune.env)
mod <- adonis(dune ~ Management + A1, data = dune.env)
## use permustats
perm <- permustats(mod)
summary(perm)
densityplot(perm)
qqmath(perm)

## example of multiple types of statistic
mod <- with(dune.env, betadisper(vegdist(dune), Management))
pmod <- permutest(mod, nperm = 99, pairwise = TRUE)
perm <- permustats(pmod)
summary(perm, interval = 0.90)
}

\keyword{ distribution }
\keyword{ smooth }
