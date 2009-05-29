\name{procrustes}
\alias{procrustes}
\alias{print.procrustes}
\alias{summary.procrustes}
\alias{print.summary.procrustes}
\alias{plot.procrustes}
\alias{points.procrustes}
\alias{text.procrustes}
\alias{lines.procrustes}
\alias{residuals.procrustes}
\alias{fitted.procrustes}
\alias{protest}
\alias{print.protest}

\title{Procrustes Rotation of Two Configurations and PROTEST }
\description{
 Function \code{procrustes} rotates a configuration to maximum similarity
  with another configuration. Function \code{protest} tests the
  non-randomness (`significance') between two configurations.
}
\usage{
procrustes(X, Y, scale = TRUE, symmetric = FALSE, scores = "sites", ...)
\method{summary}{procrustes}(object, digits = getOption("digits"), ...)
\method{plot}{procrustes}(x, kind=1, choices=c(1,2), xlab, ylab, main,
     ar.col = "blue", len=0.05, ...)
\method{points}{procrustes}(x, display = c("target", "rotated"), ...)
\method{text}{procrustes}(x, display = c("target", "rotated"), labels, ...)
\method{lines}{procrustes}(x, type = c("segments", "arrows"), choices = c(1, 2), ...)  
\method{residuals}{procrustes}(object, ...)
\method{fitted}{procrustes}(object, truemean = TRUE, ...)
protest(X, Y, scores = "sites", permutations = 999, strata, ...)
}

\arguments{
  \item{X}{Target matrix}
  \item{Y}{Matrix to be rotated.}
  \item{scale}{Allow scaling of axes of \code{Y}.}
  \item{symmetric}{Use symmetric Procrustes statistic (the rotation will
    still be non-symmetric).}
  \item{scores}{Kind of scores used. This is the \code{display} argument
  used with the corresponding \code{scores} function: see
  \code{\link{scores}}, \code{\link{scores.cca}} and
  \code{\link{scores.cca}} for alternatives.}
  \item{x, object}{An object of class \code{procrustes}.}
  \item{digits}{Number of digits in the output.}
  \item{kind}{For \code{plot} function, the kind of plot produced:
    \code{kind = 1} plots shifts in two configurations, \code{kind = 0}
    draws a corresponding empty plot, and \code{kind = 2}
    plots an impulse diagram of residuals.}
  \item{choices}{Axes (dimensions) plotted.}
  \item{xlab, ylab}{Axis labels, if defaults unacceptable.}
  \item{main}{Plot title, if default unacceptable.}
  \item{display}{Show only the \code{"target"} or \code{"rotated"}
    matrix as points.}
  \item{type}{Combine \code{target} and \code{rotated} points with line
    segments or arrows.}
  \item{truemean}{Use the original range of target matrix instead of
    centring the fitted values.}
  \item{permutations}{Number of permutation to assess the significance
    of the symmetric Procrustes statistic. }
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
  \item{ar.col}{Arrow colour.}
  \item{len}{Width of the arrow head.}
  \item{labels}{Character vector of text labels. Rownames of the result 
    object are used as default.}
  \item{...}{Other parameters passed to functions. In \code{procrustes}
    and \code{protest} parameters are passed to \code{\link{scores}}, in
    graphical functions to underlying graphical functions.}
}
\details{
  Procrustes rotation rotates a matrix to maximum similarity with a
  target matrix minimizing sum of squared differences.  Procrustes
  rotation is typically used in comparison of ordination results.  It is
  particularly useful in comparing alternative solutions in
  multidimensional scaling.  If \code{scale=FALSE}, the function only
  rotates matrix \code{Y}. If \code{scale=TRUE}, it scales linearly
  configuration \code{Y} for maximum similarity.  Since \code{Y} is scaled
  to fit \code{X}, the scaling is non-symmetric. However, with
  \code{symmetric=TRUE}, the configurations are scaled to equal
  dispersions and  a symmetric version of the Procrustes statistic
  is computed.

  Instead of matrix, \code{X} and \code{Y} can be results from an
  ordination from which \code{\link{scores}} can extract results.
  Function \code{procrustes} passes extra arguments to
  \code{\link{scores}}, \code{\link{scores.cca}} etc. so that you can
  specify arguments such as \code{scaling}. 

  Function \code{plot} plots a \code{procrustes}
  object and  returns invisibly an \code{ordiplot} object so that
  function \code{\link{identify.ordiplot}} can be used for identifying
  points. The items in the \code{ordiplot} object are called
  \code{heads} and \code{points} with \code{kind=1} (ordination diagram)
  and \code{sites} with \code{kind=2} (residuals).  In ordination
  diagrams, the arrow heads point to the target configuration, which may
  be either logical or illogical. Target and original rotated axes are
  shown as cross hairs in two-dimensional Procrustes analysis, and with
  a higher number of dimensions, the rotated axes are projected onto
  plot with their scaled and centred range. Function
  \code{plot} passes
  parameters to underlying plotting functions.  For full control of
  plots, you can draw the axes using \code{plot} with \code{kind = 0},
  and then add items with \code{points} or \code{lines}.  These
  functions pass all parameters to the underlying functions so that you
  can select the plotting characters, their size, colours etc., or you
  can select the width, colour and type of line \code{\link{segments}} or
  arrows, or you can select the orientation and head width of
  \code{\link{arrows}}.

  Function \code{residuals} returns the pointwise
  residuals, and \code{fitted} the fitted values, either centred to zero
  mean (if \code{truemean=FALSE}) or with the original scale (these
  hardly make sense if \code{symmetric = TRUE}). In
  addition, there are \code{summary} and \code{print} methods.

  If matrix \code{X} has a lower number of columns than matrix
  \code{Y}, then matrix \code{X} will be filled with zero columns to
  match dimensions. This means that the function can be used to rotate
  an ordination configuration to an environmental variable (most
  practically extracting the result with the \code{fitted} function).

  Function \code{protest} calls \code{procrustes(..., symmetric = TRUE)}
  repeatedly to estimate the `significance' of the Procrustes
  statistic. Function \code{protest} uses a correlation-like statistic
  derived from the symmetric Procrustes sum of squares \eqn{ss} as
  \eqn{r =\sqrt{(1-ss)}}, and sometimes called \eqn{m_{12}}. Function
  \code{protest} has own \code{print} method, but otherwise uses
  \code{procrustes} methods. Thus \code{plot} with a \code{protest} object
  yields a ``Procrustean superimposition plot.''
}

\value{
  Function \code{procrustes} returns an object of class
  \code{procrustes} with items. Function \code{protest} inherits from
  \code{procrustes}, but amends that with some new items:
  \item{Yrot }{Rotated matrix \code{Y}.}
  \item{X}{Target matrix.}
  \item{ss }{Sum of squared differences between \code{X} and \code{Yrot}.}
  \item{rotation}{Orthogonal rotation matrix.}
  \item{translation}{Translation of the origin.}
  \item{scale}{Scaling factor.}
  \item{symmetric}{Type of \code{ss} statistic.}
  \item{call}{Function call.}
  \item{t0}{This and the following items are only in class
    \code{protest}:  Procrustes correlation from non-permuted solution.}
  \item{t}{Procrustes correlations from permutations.}
  \item{signif}{`Significance' of \code{t}}
  \item{permutations}{Number of permutations.}
  \item{strata}{The name of the stratifying variable.}
  \item{stratum.values}{Values of the stratifying variable.}
}
\references{
  Mardia, K.V., Kent, J.T. and Bibby,
  J.M. (1979). \emph{Multivariate Analysis}. Academic Press.

  Peres-Neto, P.R. and Jackson, D.A. (2001). How well do multivariate
  data sets match? The advantages of a Procrustean superimposition
  approach over the Mantel test. \emph{Oecologia} 129: 169-178.
  
}
\author{Jari Oksanen }

\note{The function \code{protest} follows Peres-Neto & Jackson (2001),
  but the implementation is still after Mardia \emph{et al.}
  (1979).}

\seealso{\code{\link[MASS]{isoMDS}}, \code{\link{initMDS}} for obtaining
objects for \code{procrustes}, and \code{\link{mantel}} for an
alternative to \code{protest} without need of dimension reduction.} 

\examples{
data(varespec)
vare.dist <- vegdist(wisconsin(varespec))
library(MASS)  ## isoMDS
mds.null <- isoMDS(vare.dist, tol=1e-7)
mds.alt <- isoMDS(vare.dist, initMDS(vare.dist), maxit=200, tol=1e-7)
vare.proc <- procrustes(mds.alt, mds.null)
vare.proc
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)
}
\keyword{multivariate }
\keyword{htest}
