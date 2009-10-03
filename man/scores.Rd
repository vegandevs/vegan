\name{scores}
\alias{scores}
\alias{scores.default}
\alias{scores.lda}

\title{ Get Species or Site Scores from an Ordination }
\description{
  Function to access either species or site scores for specified axes in
  some ordination methods.
}
\usage{
\method{scores}{default}(x, choices, display=c("sites", "species"), ...)
}

\arguments{
  \item{x}{ An ordination result. }
  \item{choices}{ Ordination axes.  If missing, returns all axes.}
  \item{display}{ Partial match to access scores for \code{sites} or
    \code{species}.}
  \item{...}{ Other parameters (unused). }
}
\details{
  Functions \code{\link{cca}} and \code{\link{decorana}} have specific
  \code{scores} function to access their ordination scores.  Most
  standard ordination methods of libraries \pkg{mva}, \pkg{multiv} and
  \pkg{MASS} do not have a  specific \code{class}, and no specific method can be
  written for them.  However, \code{scores.default} guesses where
  some commonly used functions keep their site scores and possible
  species scores.  For site scores, the function seeks items in order
  \code{points}, \code{rproj}, \code{x}, and \code{scores}.  For species,
  the seeking order is \code{cproj}, \code{rotation}, and
  \code{loadings}.
  If \code{x} is a matrix, \code{scores.default} returns the chosen
  columns of that matrix, ignoring whether species or sites were
  requested (do not regard this as a bug but as a feature, please).
  Currently the function seems to work at least for \code{\link[MASS]{isoMDS}},
  \code{\link{prcomp}}, \code{\link{princomp}}. It may work in
  other cases or fail mysteriously.
}
\value{
  The function returns a matrix of requested scores.
}
\author{Jari Oksanen }

\seealso{\code{\link{scores.cca}}, \code{\link{scores.decorana}}.  These
have somewhat different interface -- \code{\link{scores.cca}} in
particular -- but all work with keywords \code{display="sites"} and
\code{display="species"} and return a matrix with these.
}
\examples{
data(varespec)
vare.pca <- prcomp(varespec)
scores(vare.pca, choices=c(1,2))
}
\keyword{ multivariate }




