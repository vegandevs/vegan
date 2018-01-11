\name{vegan-defunct}
\alias{vegan-defunct}

\alias{getNumObs}
\alias{permuted.index2}
\alias{metaMDSrotate}
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
\alias{commsimulator}
\alias{read.cep}

%------ NOTE:  ../R/vegan-defunct.R   must be synchronized with this!
\title{Defunct Functions in Package \pkg{vegan}}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  The functions or variables listed here are no longer part of
  \pkg{vegan} as they are no longer needed.
}
\usage{
## defunct in vegan 2.5-0
commsimulator(x, method, thin=1)

## defunct in vegan 2.4-6
read.cep(file, maxdata = 10000, positive = TRUE, trace = FALSE, force = FALSE)

## defunct in vegan 2.4-0
\method{density}{adonis}(x, ...)
\method{plot}{vegandensity}(x, main = NULL, xlab = NULL, ylab = "Density", 
   type = "l", zero.line = TRUE, obs.line = TRUE, ...)
\method{densityplot}{adonis}(x, data, xlab = "Null", ...)

## defunct in vegan 2.2-0
metaMDSrotate(object, vec, na.rm = FALSE, ...)

## defunct in vegan 2.0-0
getNumObs(object, \dots)
permuted.index2(n, control = permControl())
}

\details{

  Function \code{commsimulator} was replaced with
  \code{\link{make.commsim}} which defines the Null models, and
  functions \code{\link{nullmodel}} and
  \code{\link{simulate.nullmodel}} that check the input data and
  generate the Null model communities.

  \code{read.cep} was used to read a file formatted by relaxed strict
  CEP format used by \code{Canoco} software, among others. This
  function has been removed from \pkg{vegan}, and it is now available
  in separate package \pkg{cepreader} as function
  \code{\link[cepreader]{readCEP}}.

  The deprecated \code{density} and \code{densityplot} methods were
  replaced with similar methods for \code{\link{permustats}}. The
  \code{\link{permustats}} offers more powerful analysis tools for
  permutations, including \code{\link{summary.permustats}} giving
  \eqn{z} values (a.k.a. standardized effect sizes, SES), and Q-Q
  plots (\code{\link{qqnorm.permustats}},
  \code{\link{qqmath.permustats}}).

  Function \code{metaMDSrotate} was replaced with
  \code{\link{MDSrotate}} which can handle \code{\link{monoMDS}}
  results in addition to \code{\link{metaMDS}}.

  The permutation functions were moved to the \pkg{permute} package,
  and they are documented there.  The \pkg{permute} package replaces
  \code{permuted.index2} with \code{\link[permute]{shuffle}} and
  \code{getNumObs} with its specific
  \code{\link[permute]{nobs-methods}}.

}

\seealso{
  \code{\link{Defunct}}, \code{\link{vegan-deprecated}}
}

\keyword{internal}
