\name{vegan-defunct}
\alias{vegan-defunct}

\alias{getNumObs}
\alias{permuted.index2}
\alias{permuted.index}
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
\alias{vegandocs}
\alias{as.mlm}
\alias{as.mlm.cca}
\alias{as.mlm.rda}
\alias{humpfit}

%------ NOTE:  ../R/vegan-defunct.R   must be synchronized with this!
\title{Defunct Functions in Package \pkg{vegan}}
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  The functions or variables listed here are no longer part of
  \pkg{vegan} as they are no longer needed.
}
\usage{
## defunct in vegan 2.6-0
as.mlm(x)
humpfit(mass, spno, family = poisson, start)
vegandocs(doc = c("NEWS", "ONEWS", "FAQ-vegan", "intro-vegan",
    "diversity-vegan", "decision-vegan", "partitioning", "permutations"))

## defunct in vegan 2.5-0
commsimulator(x, method, thin=1)

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

  \code{as.mlm} function is replaced with a set functions that can
  find the same statistics directly from the ordination result object:
  see \code{\link{hatvalues.cca}}, \code{\link{rstandard.cca}},
  \code{\link{rstudent.cca}}, \code{\link{cooks.distance.cca}},
  \code{\link{vcov.cca}}.

  Function \code{humpfit} was transferred to the \pkg{natto} package and
  is still available from \url{https://github.com/jarioksa/natto/}.

  \R{} functions \code{\link{news}} should be used to read \pkg{vegan}
  NEWS (\code{news(package = "vegan")}), and
  \code{\link{browseVignettes}} is a better tool for reading vignettes
  than \code{vegandocs}.

  Function \code{commsimulator} is replaced with
  \code{\link{make.commsim}} which defines the Null models, and
  functions \code{\link{nullmodel}} and
  \code{\link{simulate.nullmodel}} that check the input data and
  generate the Null model communities.

  The deprecated \code{density} and \code{densityplot} methods are
  replaced with similar methods for \code{\link{permustats}}. The
  \code{\link{permustats}} offers more powerful analysis tools for
  permutations, including \code{\link{summary.permustats}} giving
  \eqn{z} values (a.k.a. standardized effect sizes, SES), and Q-Q
  plots (\code{\link{qqnorm.permustats}},
  \code{\link{qqmath.permustats}}).

  Function \code{metaMDSrotate} is replaced with
  \code{\link{MDSrotate}} which can handle \code{\link{monoMDS}}
  results in addition to \code{\link{metaMDS}}.

  The permutation functions were moved to the \pkg{permute} package,
  and they are documented there.  The \pkg{permute} package replaces
  \code{permuted.index} and \code{permuted.index2} with
  \code{\link[permute]{shuffle}} and \code{getNumObs} with its
  specific \code{\link[permute]{nobs-methods}}.

}

\seealso{
  \code{\link{Defunct}}, \code{\link{vegan-deprecated}}
}

\keyword{internal}
