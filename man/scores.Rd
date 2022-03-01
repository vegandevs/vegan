\name{scores}
\alias{scores}
\alias{scores.default}
\alias{scores.lda}

\title{ Get Species or Site Scores from an Ordination }

\description{
  Function to access either species or site scores for specified axes
  in some ordination methods. The \code{scores} function is generic in
  \pkg{vegan}, and \pkg{vegan} ordination functions have their own
  \code{scores} functions that are documented separately with the
  method (see e.g. \code{\link{scores.cca}},
  \code{\link{scores.metaMDS}}, \code{\link{scores.decorana}}). This
  help file documents the default \code{scores} method that is only
  used for non-\pkg{vegan} ordination objects.
}

\usage{
\method{scores}{default}(x, choices, display=c("sites", "species"), ...)
}

\arguments{
  \item{x}{ An ordination result. }
  \item{choices}{ Ordination axes.  If missing, default method returns all axes.}
  \item{display}{ Partial match to access scores for \code{sites} or
    \code{species}.}
  \item{...}{ Other parameters (unused). }
}
\details{
  Function \code{scores} is a generic method in \pkg{vegan}. Several
  \pkg{vegan} functions have their own \code{scores} methods with their
  own defaults and with some new arguments. This help page describes
  only the default method. For other methods, see, e.g.,
  \code{\link{scores.cca}}, \code{\link{scores.rda}},
  \code{\link{scores.decorana}}.

  All \pkg{vegan} ordination functions should have a \code{scores}
  method which should be used to extract the scores instead of
  directly accessing them. Scaling and transformation of scores should
  also happen in the \code{scores} function. If the \code{scores}
  function is available, the results can be plotted using
  \code{\link{ordiplot}}, \code{\link{ordixyplot}} etc., and the
  ordination results can be compared in \code{\link{procrustes}}
  analysis.
  
  The \code{scores.default} function is used to extract scores from
  non-\pkg{vegan} ordination results.  Many standard ordination
  methods of libraries do not have a specific \code{class}, and no
  specific method can be written for them.  However,
  \code{scores.default} guesses where some commonly used functions
  keep their site scores and possible species scores.

  If \code{x} is a matrix, \code{scores.default} returns the chosen
  columns of that matrix, ignoring whether species or sites were
  requested (do not regard this as a bug but as a feature, please).
  Currently the function seems to work at least for \code{\link[MASS]{isoMDS}},
  \code{\link{prcomp}}, \code{\link{princomp}} and some \pkg{ade4} objects. 
  It may work in other cases or fail mysteriously.
}
\value{
  The function returns a matrix of scores.
}
\author{Jari Oksanen }

\seealso{
  Specific \code{scores} functions include (but are not limited to)
  \code{\link{scores.cca}}, \code{\link{scores.rda}},
  \code{\link{scores.decorana}}, \code{\link{scores.envfit}},
  \code{\link{scores.metaMDS}}, \code{\link{scores.monoMDS}} and
  \code{\link{scores.pcnm}}.  These have somewhat different interface
  -- \code{\link{scores.cca}} in particular -- but all work with
  keywords \code{display="sites"} and return a matrix. However, they
  may also return a list of matrices, and some other \code{scores}
  methods will have quite different arguments.
}

\examples{
data(varespec)
vare.pca <- prcomp(varespec)
scores(vare.pca, choices=c(1,2))
}
\keyword{ multivariate }




