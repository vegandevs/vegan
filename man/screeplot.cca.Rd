\name{screeplot.cca}
\alias{screeplot.cca}
\alias{screeplot.princomp}
\alias{screeplot.prcomp}
\alias{screeplot.decorana}
\alias{bstick}
\alias{bstick.default}
\alias{bstick.cca}
\alias{bstick.prcomp}
\alias{bstick.princomp}
\alias{bstick.decorana}
\title{Screeplots for Ordination Results and Broken Stick Distributions}
\description{
  Screeplot methods for plotting variances of ordination axes/components
  and overlaying broken stick distributions. Also, provides alternative
  screeplot methods for \code{\link{princomp}} and \code{\link{prcomp}}.
}
\usage{
\method{screeplot}{cca}(x, bstick = FALSE, type = c("barplot", "lines"),
         npcs = min(10, if (is.null(x$CCA) || x$CCA$rank == 0) x$CA$rank else x$CCA$rank),
         ptype = "o", bst.col = "red", bst.lty = "solid",
         xlab = "Component", ylab = "Inertia",
         main = deparse(substitute(x)), legend = bstick,
         \dots)

\method{screeplot}{decorana}(x, bstick = FALSE, type = c("barplot", "lines"),
         npcs = 4,
         ptype = "o", bst.col = "red", bst.lty = "solid",
         xlab = "Component", ylab = "Inertia",
         main = deparse(substitute(x)), legend = bstick,
         \dots)

\method{screeplot}{prcomp}(x, bstick = FALSE, type = c("barplot", "lines"),
         npcs = min(10, length(x$sdev)),
         ptype = "o", bst.col = "red", bst.lty = "solid",
         xlab = "Component", ylab = "Inertia",
         main = deparse(substitute(x)), legend = bstick,
         \dots)

\method{screeplot}{princomp}(x, bstick = FALSE, type = c("barplot", "lines"),
         npcs = min(10, length(x$sdev)),
         ptype = "o", bst.col = "red", bst.lty = "solid",
         xlab = "Component", ylab = "Inertia",
         main = deparse(substitute(x)), legend = bstick,
         \dots)

bstick(n, \dots)

\method{bstick}{default}(n, tot.var = 1, \dots)

\method{bstick}{cca}(n, \dots)

\method{bstick}{prcomp}(n, \dots)

\method{bstick}{princomp}(n, \dots)

\method{bstick}{decorana}(n, \dots)
}
\arguments{
  \item{x}{an object from which the component variances can be determined.}
  \item{bstick}{logical; should the broken stick distribution be drawn?}
  \item{npcs}{the number of components to be plotted.}
  \item{type}{the type of plot.}
  \item{ptype}{if \code{type == "lines"} or \code{bstick = TRUE}, a
    character indicating the type of plotting used for the lines;
    actually any of the \code{type}s as in \code{\link{plot.default}}.} 
  \item{bst.col, bst.lty}{the colour and line type used to draw the
    broken stick distribution.}
  \item{xlab, ylab, main}{graphics parameters.}
  \item{legend}{logical; draw a legend?}
  \item{n}{an object from which the variances can be extracted or the
    number of variances (components) in the case of
    \code{bstick.default}.}
  \item{tot.var}{the total variance to be split.}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  The functions provide screeplots for most ordination methods in
  \pkg{vegan} and enhanced versions with broken stick for
  \code{\link{prcomp}} and \code{\link{princomp}}. 

  Function \code{bstick} gives the brokenstick values which are ordered
  random proportions, defined as  \eqn{p_i = (tot/n) \sum_{x=i}^n 
    (1/x)}{p[i] = tot/n sum(from x=i to n) 1/x} (Legendre & Legendre 2012), where
  \eqn{tot} is the total  and \eqn{n} is the number of brokenstick
  components (cf. \code{\link{radfit}}).  Broken stick has
  been recommended as a stopping rule in principal component analysis
  (Jackson 1993): principal components should be retained as long as
  observed eigenvalues are higher than corresponding random broken stick
  components.

  The \code{bstick} function is generic. The default needs the number of
  components and the total, and specific methods extract this
  information from ordination results.  There also is a \code{bstick}
  method for \code{\link{cca}}.  However, the broken stick model is not
  strictly valid for correspondence analysis (CA), because eigenvalues
  of CA are defined to be \eqn{\leq 1}{<=1}, whereas brokenstick
  components have no such restrictions. The brokenstick components in
  detrended correspondence analysis (DCA) assume that input data are of
  full rank, and additive eigenvalues are used in \code{screeplot} (see
  \code{\link{decorana}}).
}

\value{
  Function \code{screeplot} draws a plot on the currently active device, 
  and returns invisibly the \code{\link{xy.coords}} of the points or
  bars for the eigenvalues.

  Function \code{bstick} returns a numeric vector of broken stick
  components. 
}
\references{
  Jackson, D. A. (1993). Stopping rules in principal components
  analysis: a comparison of heuristical and statistical
  approaches. \emph{Ecology} 74, 2204--2214.

  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd English
  ed. Elsevier.
  }

\author{Gavin L. Simpson}
\seealso{
  \code{\link{cca}}, \code{\link{decorana}}, \code{\link{princomp}} and
  \code{\link{prcomp}} for the ordination functions, and
  \code{\link{screeplot}} for the stock version.  }
\examples{
data(varespec)
vare.pca <- rda(varespec, scale = TRUE)
bstick(vare.pca)
screeplot(vare.pca, bstick = TRUE, type = "lines")
}
\keyword{multivariate}
