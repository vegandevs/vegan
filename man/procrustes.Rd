\name{procrustes}
\alias{procrustes}
\alias{summary.procrustes}
\alias{plot.procrustes}
\alias{points.procrustes}
\alias{text.procrustes}
\alias{lines.procrustes}
\alias{residuals.procrustes}
\alias{fitted.procrustes}
\alias{predict.procrustes}
\alias{protest}

\title{Procrustes Rotation of Two Configurations and PROTEST }
\description{
 Function \code{procrustes} rotates a configuration to maximum similarity
  with another configuration. Function \code{protest} tests the
  non-randomness (`significance') between two configurations.
}
\usage{
procrustes(X, Y, scale = TRUE, symmetric = FALSE, scores = "sites", ...)
\method{summary}{procrustes}(object, digits = getOption("digits"), ...)
\method{plot}{procrustes}(x, kind=1, choices=c(1,2), to.target = TRUE, 
    type = "p", xlab, ylab, main, ar.col = "blue", len=0.05, 
    cex = 0.7, ...)
\method{points}{procrustes}(x, display = c("target", "rotated"), ...)
\method{text}{procrustes}(x, display = c("target", "rotated"), labels, ...)
\method{lines}{procrustes}(x, type = c("segments", "arrows"), choices = c(1, 2), ...)  
\method{residuals}{procrustes}(object, ...)
\method{fitted}{procrustes}(object, truemean = TRUE, ...)
\method{predict}{procrustes}(object, newdata, truemean = TRUE, ...)
protest(X, Y, scores = "sites", permutations = how(nperm = 999), ...)
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
  \item{to.target}{Draw arrows to point to target.}
  \item{type}{The type of plot drawn. In \code{plot}, the \code{type}
    can be \code{"points"} or \code{"text"} to select the marker for
    the tail of the arrow, or \code{"none"} for drawing an empty
    plot. In \code{lines} the \code{type} selects either
    \code{\link{arrows}} or line \code{\link{segments}} to connect
    target and rotated configuration.}
  \item{truemean}{Use the original range of target matrix instead of
    centring the fitted values. Function \code{plot.procrustes} needs
    \code{truemean = FALSE}.}
  \item{newdata}{Matrix of coordinates to be rotated and translated to
     the target.}
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{ar.col}{Arrow colour.}
  \item{len}{Width of the arrow head.}
  \item{labels}{Character vector of text labels. Rownames of the result 
    object are used as default.}
  \item{cex}{Character expansion for points or text.}
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

  Function \code{plot} plots a \code{procrustes} object and returns
  invisibly an \code{ordiplot} object so that function
  \code{\link{identify.ordiplot}} can be used for identifying
  points. The items in the \code{ordiplot} object are called
  \code{heads} and \code{points} with \code{kind=1} (ordination
  diagram) and \code{sites} with \code{kind=2} (residuals).  In
  ordination diagrams, the arrow heads point to the target
  configuration if \code{to.target = TRUE}, and to rotated
  configuration if \code{to.target = FALSE}.  Target and original
  rotated axes are shown as cross hairs in two-dimensional Procrustes
  analysis, and with a higher number of dimensions, the rotated axes
  are projected onto plot with their scaled and centred
  range. Function \code{plot} passes parameters to underlying plotting
  functions.  For full control of plots, you can draw the axes using
  \code{plot} with \code{kind = 0}, and then add items with
  \code{points} or \code{lines}.  These functions pass all parameters
  to the underlying functions so that you can select the plotting
  characters, their size, colours etc., or you can select the width,
  colour and type of line \code{\link{segments}} or arrows, or you can
  select the orientation and head width of \code{\link{arrows}}.

  Function \code{residuals} returns the pointwise
  residuals, and \code{fitted} the fitted values, either centred to zero
  mean (if \code{truemean=FALSE}) or with the original scale (these
  hardly make sense if \code{symmetric = TRUE}). In
  addition, there are \code{summary} and \code{print} methods.

  If matrix \code{X} has a lower number of columns than matrix
  \code{Y}, then matrix \code{X} will be filled with zero columns to
  match dimensions. This means that the function can be used to rotate
  an ordination configuration to an environmental variable (most
  practically extracting the result with the \code{fitted}
  function). Function \code{predict} can be used to add new rotated
  coordinates to the target. The \code{predict} function will always
  translate coordinates to the original non-centred matrix. The
  function cannot be used with \code{newdata} for \code{symmetric}
  analysis.

  Function \code{protest} performs symmetric Procrustes analysis
  repeatedly to estimate the `significance' of the Procrustes
  statistic. Function \code{protest} uses a correlation-like statistic
  derived from the symmetric Procrustes sum of squares \eqn{ss} as
  \eqn{r =\sqrt{1-ss}}{r = sqrt(1-ss)}, and also prints the sum of
  squares of the symmetric analysis, sometimes called
  \eqn{m_{12}^2}{squared m12}. Function \code{protest} has own
  \code{print} method, but otherwise uses \code{procrustes}
  methods. Thus \code{plot} with a \code{protest} object yields a
  ``Procrustean superimposition plot.''  }

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
  \item{xmean}{The centroid of the target.}
  \item{symmetric}{Type of \code{ss} statistic.}
  \item{call}{Function call.}
  \item{t0}{This and the following items are only in class
    \code{protest}:  Procrustes correlation from non-permuted solution.}
  \item{t}{Procrustes correlations from permutations. The distribution
    of these correlations can be inspected with \code{\link{permustats}}
    function.}
  \item{signif}{`Significance' of \code{t}}
  \item{permutations}{Number of permutations.}
  \item{control}{A list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}.}
  \item{control}{the list passed to argument \code{control} describing
    the permutation design.}
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

\seealso{\code{\link{monoMDS}},  for obtaining
objects for \code{procrustes}, and \code{\link{mantel}} for an
alternative to \code{protest} without need of dimension reduction. See
\code{\link[permute]{how}} for details on specifying the type of
permutation required.} 

\examples{
data(varespec)
vare.dist <- vegdist(wisconsin(varespec))
mds.null <- monoMDS(vare.dist, y = cmdscale(vare.dist))
mds.alt <- monoMDS(vare.dist)
vare.proc <- procrustes(mds.alt, mds.null)
vare.proc
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)
}
\keyword{multivariate }
\keyword{htest}
