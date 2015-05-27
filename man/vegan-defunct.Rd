\name{vegan-defunct}

\alias{density.adonis}
\alias{density.anosim}
\alias{density.mantel}
\alias{density.mrpp}
\alias{density.permutest.cca}
\alias{density.protest}
\alias{plot.vegandensity}
\alias{densityplot.adonis}
\alias{density.oecosimu}
\alias{densityplot.oecosimu}
\alias{vegan-defunct}
%------ NOTE:  ../R/vegan-defunct.R   must be synchronized with this!
\title{Defunct Functions in Package \pkg{vegan}}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  The functions or variables listed here are no longer part of
  \pkg{vegan} as they are no longer needed.
}
\usage{
\method{density}{adonis}(x, ...)
\method{plot}{vegandensity}(x, main = NULL, xlab = NULL, ylab = "Density", 
   type = "l", zero.line = TRUE, obs.line = TRUE, ...)
\method{densityplot}{adonis}(x, data, xlab = "Null", ...)
}
%\arguments{
%%  \item{x}{An object with permutation results.}
%  \item{main, xlab, ylab, type, zero.line}{Arguments of
%    \code{\link{plot.density}}, \code{\link[lattice]{densityplot}}.}
%  \item{obs.line}{Draw vertical line for the observed
%    statistic. Logical value \code{TRUE} draws a red line, and
%    \code{FALSE} draws nothing. Alternatively, \code{obs.line} can be a
%    definition of the colour used for the line, either as a numerical
%%    value from the \code{\link[grDevices]{palette}} or as the name of
%    the colour, or other normal definition of the colour.}
%  \item{data}{Ignored.}
%}


\details{
  The deprecated \code{density} and \code{densityplot} methods are
  replaced with similar methods for \code{\link{permustats}}. The
  \code{\link{permustats}} offers more powerful analysis tools for
  permutations, including \code{\link{summary.permustats}} giving
  \eqn{z} values (a.k.a. standardized effect sizes, SES), and Q-Q
  plots (\code{\link{qqnorm.permustats}},
  \code{\link{qqmath.permustats}}.
}

\seealso{
  \code{\link{Defunct}}, \code{\link{vegan-deprecated}}
}

\keyword{internal}
